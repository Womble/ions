(load "/home/tomd/programs/lisp/maths.lisp")
(proclaim '(inline scalar-pot vector-pot rot-mat q-add))
(proclaim '(optimize (speed 3) (compilation-speed 0) (debug 3)))



(defconstant epsi0 8.85418782d-12)
(defconstant mu0  (/ 1 c0 c0 epsi0))
(defconstant e-charge -1.60217646d-19)
(defconstant p-mass 1.67262158d-27)
(defconstant e-mass 9.10938188d-31)
(defconstant delta-t 1d-3)
(defconstant length-scale 5d-1)
(defconstant +n-ions+ 20)

(defparameter *ion-array*
  (make-array `(,+n-ions+ 8) :element-type 'double-float))
(defparameter *grid* 
  (make-array '(7 7 7) :initial-element '(0 0 0 0 0 0 0 0 0) :element-type 'list))

(defmacro posx (n)   `(the double-float (aref *ion-array* ,n 0)))
(defmacro posy (n)   `(the double-float (aref *ion-array* ,n 1)))
(defmacro posz (n)   `(the double-float (aref *ion-array* ,n 2)))
(defmacro velx (n)   `(the double-float (aref *ion-array* ,n 3)))
(defmacro vely (n)   `(the double-float (aref *ion-array* ,n 4)))
(defmacro velz (n)   `(the double-float (aref *ion-array* ,n 5)))
(defmacro charge (n) `(the double-float (aref *ion-array* ,n 6)))
(defmacro mass (n)   `(the double-float (aref *ion-array* ,n 7)))

(dotimes (i +n-ions+)
  (let ((r (coerce (gauss-random length-scale (* length-scale 3)) 'double-float))
	(theta (coerce (random (* 2 pi)) 'double-float))
	(z (coerce (gauss-random length-scale) 'double-float))
	(charge (* 1d0 (expt -1d0 (mod i 2)))))
    (setf (posx i) (coerce (* r (cos theta)) 'double-float)
	  (posy i) (coerce (* r (sin theta)) 'double-float)
	  (posz i) z
	  (velx i) (* 1d-5 (- 1d0 (random 3)))
	  (vely i) (* 1d-5 (- 1d0 (random 3)))
	  (velz i) (* 1d-5 (- 1d0 (random 3)))
	  (charge i) charge
	  (mass i) (if (= charge -1) 1d0 1839d0))))


(defun scalar-pot (x y z &optional (self -1))
  (declare (type double-float x y z)
	   (type fixnum self))
  (let ((ans 0d0))
    (dotimes (i +n-ions+)
      (if (= i self)()
	  (decf ans (/ (charge i)
		       (the double-float (q-add (- x (posx i)) (- y (posy i)) (- z (posz i))))))))
    	ans))

(defun vector-pot (x y z &optional (self -1))
  (declare (type double-float x y z)
	   (type fixnum self))
  (let ((ax 0d0)(ay 0d0)(az 0d0))
    (dotimes (i +n-ions+)
      (if (= i self) ()
	  (let ((r (the double-float (q-add (- x (posx i)) (- y (posy i)) (- z (posz i))))))
	    (incf ax (/ (* (charge i) (velx i) ) r))
	    (incf ay (/ (* (charge i) (vely i) ) r))
	    (incf az (/ (* (charge i) (velx i) ) r)))))
	(list ax ay az)))

(defun s-grad (fn)
  "gradient of 3D scalar field given by FN"
  #'(lambda (x y z &optional (self -1))
      (declare (type double-float x y z))
      (list (funcall (partial-diff fn 0 4) x y z self)
	    (funcall (partial-diff fn 1 4) x y z self)
	    (funcall (partial-diff fn 2 4) x y z self))))

(defun s-curl (fn)
  "curl of 3D vector field given by FN"
  (let ((d/dx (partial-diff-vec fn 0 4))
	(d/dy (partial-diff-vec fn 1 4))
	(d/dz (partial-diff-vec fn 2 4)))
  #'(lambda (x y z &optional (self -1))
      (declare (type double-float x y z))
      (let ((ddx (funcall d/dx x y z self))
	    (ddy (funcall d/dy x y z self))
	    (ddz (funcall d/dz x y z self)))
	(list (- (the double-float (caddr ddy))
		 (the double-float (cadr ddz)))
	      (- (the double-float (car ddz))
		 (the double-float (caddr ddx)))
	      (- (the double-float (cadr ddx))
		 (the double-float (car ddy))))))))

(defun q-add (x y z)
  (declare (type double-float x y z))
  (sqrt (+ (* x x) (* y y) (* z z))))

(defun update-grid (&optional (ext-B-field #'(lambda (x y z) (declare (ignore x y z)) '(0 0 0))))
  (dotimes (i 7)
    (dotimes (j 7)
      (dotimes (k 7)
	(let* ((pos (list (coerce (- (* i 2/3) 2) 'double-float) (coerce (- (* j 2/3) 2) 'double-float) (coerce (- (* k 2/3) 2) 'double-float)))
	       (E-field (s-grad (memoize #'scalar-pot)))
	       (B-field (s-curl (memoize #'vector-pot)))
	       (B (vec+ (apply ext-b-field pos) (apply B-field pos)))
	       (E (apply E-field pos)))
	  (setf (aref *grid* i j k) (append E B Pos)))))))

(defun interpolate (pos)
  (dbind (i j k) (mapcar #[* 3/2 (+ 2 $)] pos)
    (multiple-value-bind (i_a i_b) (floor i)    ;i_a is the x reference of the cell in position below the position (i.e the first cell with a position < the position) 
      (multiple-value-bind (j_a j_b) (floor j)  ;i_b is the fraction of the way it it between the inferior grid reference and the superior one
	(multiple-value-bind (k_a k_b) (floor k)
	  (if (or (< i 0) (>= i 6) (< j 0) (>= j 6) (< k 0) (>= k 6))
	      '(0d0 0d0 0d0 0d0 0d0 0d0)
	      (dbind ((exa eya eza bxa bya bza) (exb eyb ezb bxb byb bzb))
		  (list (aref *grid* i_a j_a k_a) (aref *grid* (1+ i_a) (1+ j_a) (1+ k_a)))
		(list (+ (* exa i_b) (* exb (- 1 i_b))) (+ (* eya j_b) (* eyb (- 1 j_b))) (+ (* eza k_b) (* ezb (- 1 k_b)))
		      (+ (* bxa i_b) (* bxb (- 1 i_b))) (+ (* bya j_b) (* byb (- 1 j_b))) (+ (* bza k_b) (* bzb (- 1 k_b)))))))))))

  

;;;readouts
(defun positions (stream)
  (dotimes (i +n-ions+)
    (format stream "~&~s ~s ~s~&" (coerce (posx i) 'single-float)(coerce (posy i) 'single-float)(coerce (posz i) 'single-float))))

(defun velocities (stream)
  (dotimes (i +n-ions+)
    (format stream "~&~s ~s ~s~&" (coerce (velx i) 'single-float)(coerce (vely i) 'single-float)(coerce (velz i) 'single-float))))

(defun x-v-vx (stream)
  (dotimes (i +n-ions+)
    (format stream "~&~s ~s~&" (coerce (posx i) 'single-float)(coerce (velx i) 'single-float))))

(defun E-electro (stream)
  (let ((ans 0))
    (dotimes (i +n-ions+)
      (incf ans (* (Charge i) (scalar-pot (posx i)(posy i)(posz i)i) 0.5)))
    (format stream "~S~&" (coerce ans 'single-float))))

(defun E-kinetic (stream)
  (let ((ans 0))
    (dotimes (i +n-ions+)
      (incf ans (* (mass i) (expt (q-add (velx i)(vely i)(velz i)) 2) 0.5)))
    (format stream "~S~&" (coerce ans 'single-float))))



(defun main (&optional (ticks 100))
  (with-open-file (stream "dat/dat-start-pos" :direction :output :if-exists :supersede :if-does-not-exist :create)
    (positions stream))
  (with-open-file (stream "dat/dat-start-vel" :direction :output :if-exists :supersede :if-does-not-exist :create)
    (velocities stream))
  (with-open-file (stream "dat/dat-start-x-vx" :direction :output :if-exists :supersede :if-does-not-exist :create)
    (x-v-vx stream))

  (with-open-file (stream1 "dat/dat1"   :direction :output :if-exists :supersede :if-does-not-exist :create)
    (with-open-file (stream2 "dat/dat2"   :direction :output :if-exists :supersede :if-does-not-exist :create)
      (with-open-file (stream3 "dat/dat3"   :direction :output :if-exists :supersede :if-does-not-exist :create)
	(with-open-file (stream4 "dat/dat4"   :direction :output :if-exists :supersede :if-does-not-exist :create)
	  (with-open-file (streamEE "dat/datEEs" :direction :output :if-exists :supersede :if-does-not-exist :create)
	    (with-open-file (streamEK "dat/datEK"  :direction :output :if-exists :supersede :if-does-not-exist :create)

	      
	      (dotimes (i ticks)
		
		(format streamEK "~&~S " i)(E-kinetic streamEK)
		(format streamEE "~&~S " i)(E-electro streamEE)
		(format stream4 "~&~S ~S ~S~&" (coerce (posx 4) 'single-float)(coerce (posy 4) 'single-float)(coerce (posz 4)'single-float))
		(format stream3 "~&~S ~S ~S~&" (coerce (posx 3) 'single-float)(coerce (posy 3) 'single-float)(coerce (posz 3)'single-float))
		(format stream2 "~&~S ~S ~S~&" (coerce (posx 2) 'single-float)(coerce (posy 2) 'single-float)(coerce (posz 2)'single-float))
		(format stream1 "~&~S ~S ~S~&" (coerce (posx 1) 'single-float)(coerce (posy 1) 'single-float)(coerce (posz 1)'single-float))
    
    
		(update-grid)
		(dotimes (n +n-ions+)
		  (incf (posx n) (* (velx n) delta-t))
		  (incf (posy n) (* (vely n) delta-t))
		  (incf (posz n) (* (velz n) delta-t))
		  (let ((pos (list (posx n)(posy n)(Posz n))))
		    (dbind (ex ey ez bx by bz) (interpolate pos)
		      (incf (velx n) (* 0.5d0 (the double-float ex) delta-t (/ (charge n) (mass n))))
		      (incf (vely n) (* 0.5d0 (the double-float ey) delta-t (/ (charge n) (mass n)))) ;half accel from E
		      (incf (velz n) (* 0.5d0 (the double-float ez) delta-t (/ (charge n) (mass n))))
		      (dbind ((x y z)) (matrix-multiply (cons (list (velx n) (vely n) (velz n)) nil) 
							(rot-mat (list Bx by bz) (* delta-t (/ (charge n) (mass n)) (q-add Bx by bz))))
			(setf (velx n) x
			      (vely n) y   ;rotation from B field
			      (velz n) z))
		      (incf (velx n) (* 0.5d0 (the double-float ex) delta-t (/ (charge n) (mass n))))
		      (incf (vely n) (* 0.5d0 (the double-float ey) delta-t (/ (charge n) (mass n)))) ;half accel from E
		      (incf (velz n) (* 0.5d0 (the double-float ez) delta-t (/ (charge n) (mass n))))))))))))))

  (with-open-file (stream "dat/dat-final-pos" :direction :output :if-exists :supersede :if-does-not-exist :create)
    (positions stream))
  (with-open-file (stream "dat/dat-final-vel" :direction :output :if-exists :supersede :if-does-not-exist :create)
    (velocities stream))
  (with-open-file (stream "dat/dat-final-x-vx" :direction :output :if-exists :supersede :if-does-not-exist :create)
    (x-v-vx stream)))
