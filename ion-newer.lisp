(load "/home/tomd/programs/lisp/maths.lisp")
(proclaim '(inline scalar-pot vector-pot rot-mat q-add))
(proclaim '(optimize (speed 3) (compilation-speed 0) (debug 3)))

;;low ion but no field interpolation version

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
	(charge (* 1d0 (expt -1d0 (random 2)))))
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
  #'(lambda (x y z self)
      (list (funcall (partial-diff fn 0 4) x y z self)
	    (funcall (partial-diff fn 1 4) x y z self)
	    (funcall (partial-diff fn 2 4) x y z self))))

(defun s-curl (fn)
  "curl of 3D vector field given by FN"
  (let ((d/dx (partial-diff-vec fn 0 4))
	(d/dy (partial-diff-vec fn 1 4))
	(d/dz (partial-diff-vec fn 2 4)))
  #'(lambda (x y z self)
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

(defun rot-mat (vec theta)
  "creates the rotation matrix for rotation of angle theta about vec"
  (let* ((cosin (coerce (cos theta) 'double-float))
	 (sine  (coerce (sin theta) 'double-float))
	 (1-cos (coerce (- 1 cosin) 'double-float))
	 (vec-normalised (vec-norm vec))
	 (ux (coerce (car   vec-normalised) 'double-float))
	 (uy (coerce (cadr  vec-normalised) 'double-float))
	 (uz (coerce (caddr vec-normalised) 'double-float)))
    `((   ,(+ cosin (* ux ux 1-cos))    ,(- (* ux uy 1-cos) (* uz sine))  ,(+ (* ux uz 1-cos) (* uy sine) ))
      (,(+ (* uy ux 1-cos) (* uz sine))    ,(+ cosin (* uy uy 1-cos))     ,(- (* uy uz 1-cos) (* ux sine) ))
      (,(- (* uz ux 1-cos) (* uy sine)) ,(+ (* uz uy 1-cos) (* ux sine))     ,(+ cosin (* uz uz 1-cos)    ))) ))

(defun main (&optional (ticks 100))
  (with-open-file (stream "dat-start-pos" :direction :output :if-exists :supersede :if-does-not-exist :create)
    (positions stream))
  (with-open-file (stream "dat-start-vel" :direction :output :if-exists :supersede :if-does-not-exist :create)
    (velocities stream))
  (with-open-file (stream "dat-start-x-vx" :direction :output :if-exists :supersede :if-does-not-exist :create)
    (x-v-vx stream))

  (dotimes (i ticks)

    (with-open-file (stream "dat1" :direction :output :if-exists :append :if-does-not-exist :create)
      (format stream "~&~S ~S ~S~&" (coerce (posx 1) 'single-float)(coerce (posy 1) 'single-float)(coerce (posz 1)'single-float)))
    (with-open-file (stream "dat2" :direction :output :if-exists :append :if-does-not-exist :create)
      (format stream "~&~S ~S ~S~&" (coerce (posx 2) 'single-float)(coerce (posy 2) 'single-float)(coerce (posz 2)'single-float)))
    (with-open-file (stream "dat3" :direction :output :if-exists :append :if-does-not-exist :create)
      (format stream "~&~S ~S ~S~&" (coerce (posx 3) 'single-float)(coerce (posy 3) 'single-float)(coerce (posz 3)'single-float)))
    (with-open-file (stream "dat4" :direction :output :if-exists :append :if-does-not-exist :create)
      (format stream "~&~S ~S ~S~&" (coerce (posx 4) 'single-float)(coerce (posy 4) 'single-float)(coerce (posz 4)'single-float)))
    (with-open-file (stream "datEEs" :direction :output :if-exists :append :if-does-not-exist :create)
      (format stream "~&~S " i)(E-electro stream))
    (with-open-file (stream "datEK" :direction :output :if-exists :append :if-does-not-exist :create)
      (format stream "~&~S " i)(E-kinetic stream))

    (dotimes (n +n-ions+)
      (incf (posx n) (* (velx n) delta-t))
      (incf (posy n) (* (vely n) delta-t))
      (incf (posz n) (* (velz n) delta-t))
      (let ((pos (list (posx n)(posy n)(Posz n) n))
	    (E-field (s-grad (memoize #'scalar-pot)))
	    (B-field (s-curl (memoize #'vector-pot))))
	(let ((B (vec+ '(0 0 1) (apply B-field pos)))
	      (E (apply E-field pos)))
	  (incf (velx n) (* 0.5d0 (the double-float (nth 0 E)) delta-t (/ (charge n) (mass n))))
	  (incf (vely n) (* 0.5d0 (the double-float (nth 1 E)) delta-t (/ (charge n) (mass n)))) ;half accel from E
	  (incf (velz n) (* 0.5d0 (the double-float (nth 2 E)) delta-t (/ (charge n) (mass n))))
	  (dbind ((x y z)) (matrix-multiply (cons (list (velx n) (vely n) (velz n)) nil) (rot-mat B (* delta-t (/ (charge n) (mass n)) (apply #'q-add B))))
	    (setf (velx n) x
		  (vely n) y
		  (velz n) z))
	  (incf (velx n) (* 0.5d0 (the double-float (nth 0 E)) delta-t (/ (charge n) (mass n))))
	  (incf (vely n) (* 0.5d0 (the double-float (nth 1 E)) delta-t (/ (charge n) (mass n)))) ;half accel from E
	  (incf (velz n) (* 0.5d0 (the double-float (nth 2 E)) delta-t (/ (charge n) (mass n)))))))))


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
