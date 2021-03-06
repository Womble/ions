(load "/home/tomd/programs/lisp/maths.lisp")
(proclaim '(inline scalar-pot vector-pot rot-mat q-add e-field b-field))
(proclaim '(optimize (speed 3) (compilation-speed 0) (debug 0)))


(defconstant epsi0 8.85418782d-12)
(defconstant mu0  (/ 1 c0 c0 epsi0))
(defconstant e-charge -1.60217646d-19)
(defconstant p-mass 1.67262158d-27)
(defconstant e-mass 9.10938188d-31)
(defconstant delta-t 1d-3)
(defconstant length-scale 5d-1)
(defconstant +n-ions+ (the fixnum 1000))
(defconstant +grid-points+ (the fixnum 7))

(defparameter *ion-array*
  (make-array `(,+n-ions+ 8) :element-type 'double-float :initial-element 0d0))
(defparameter *grid* 
  (make-array `(,+grid-points+ ,+grid-points+ ,+grid-points+) 
;	      :initial-element '(0d0 0d0 0d0   ;E 
;	        		 0d0 0d0 0d0   ;B
;	        		 0d0 0d0 0d0   ;Pos
;	        		 0d0           ;Rho
;	        		 0d0 0d0 0d0)  ;J	      
	      :element-type 'list))

(defmacro posx (n)   `(the double-float (aref *ion-array* ,n 0)))
(defmacro posy (n)   `(the double-float (aref *ion-array* ,n 1)))
(defmacro posz (n)   `(the double-float (aref *ion-array* ,n 2)))
(defmacro velx (n)   `(the double-float (aref *ion-array* ,n 3)))
(defmacro vely (n)   `(the double-float (aref *ion-array* ,n 4)))
(defmacro velz (n)   `(the double-float (aref *ion-array* ,n 5)))
(defmacro charge (n) `(the double-float (aref *ion-array* ,n 6)))
(defmacro mass (n)   `(the double-float (aref *ion-array* ,n 7)))


(defmacro g-Ex  (i j k) `(the double-float (nth 0 (aref *grid* ,i ,j ,k))))
(defmacro g-Ey  (i j k) `(the double-float (nth 1 (aref *grid* ,i ,j ,k))))
(defmacro g-Ez  (i j k) `(the double-float (nth 2 (aref *grid* ,i ,j ,k))))
(defmacro g-Bx  (i j k) `(the double-float (nth 3 (aref *grid* ,i ,j ,k))))
(defmacro g-By  (i j k) `(the double-float (nth 4 (aref *grid* ,i ,j ,k))))
(defmacro g-Bz  (i j k) `(the double-float (nth 5 (aref *grid* ,i ,j ,k))))
(defmemo g-pos  (i j k)  (firstn 3 (nthcdr 6 (aref *grid* i j k))))
(defmacro g-rho (i j k) `(the double-float (nth 9  (aref *grid* ,i ,j ,k))))
(defmacro g-J   (i j k) `(nthcdr 10 (aref *grid* ,i ,j ,k)))
(defmacro g-Jx  (i j k) `(the double-float (nth 10 (aref *grid* ,i ,j ,k))))
(defmacro g-Jy  (i j k) `(the double-float (nth 11 (aref *grid* ,i ,j ,k))))
(defmacro g-Jz  (i j k) `(the double-float (nth 12 (aref *grid* ,i ,j ,k))))

;(dotimes (i +n-ions+)
;  (let ((r (coerce (gauss-random length-scale (* length-scale 3)) 'double-float))
;	(theta (coerce (random (* 2 pi)) 'double-float))
;	(z (coerce (gauss-random length-scale) 'double-float))
;	(charge (* 1d0 (expt -1d0 (mod i 2)))))
;    (setf (posx i) (coerce (* r (cos theta)) 'double-float)
;	  (posy i) (coerce (* r (sin theta)) 'double-float)
;	  (posz i) z
;	  (velx i) (/ (- 1d0 (random 3d0)) 3)
;	  (vely i) (/ (- 1d0 (random 3d0)) 3)
;	  (velz i) (/ (- 1d0 (random 3d0)) 3)
;	  (charge i) charge
;	  (mass i) (if (= charge -1) 1d0 1839d0))))

;setup ion inital conditions
(dotimes (i +n-ions+)
  (let ((x (coerce (gauss-random (/ length-scale 3)) 'double-float))
	(y (coerce (gauss-random (/ length-scale 3)) 'double-float))
	(z (coerce (gauss-random (/ length-scale 3)) 'double-float))
	(charge (* 1d0 (expt -1d0 (mod i 2)))))
    (setf (posx i) x
	  (posy i) y
	  (posz i) (- z 1)
	  (velx i) (coerce (gauss-random (/ length-scale 3)) 'double-float)
	  (vely i) (coerce (gauss-random (/ length-scale 3)) 'double-float)
	  (velz i) (+ 2 (coerce (gauss-random (/ length-scale 3)) 'double-float))
	  (charge i) charge
	  (mass i) (if (= charge -1) 1d0 1839d0))))

;set up grid co-ords
(dotimes (i +grid-points+)
    (dotimes (j +grid-points+)
      (dotimes (k +grid-points+)
	(setf (aref *grid* i j k) (list 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0))
	(setf (nth 6 (aref *grid* i j k)) (coerce (- (* i (/ (* 3 length-scale))) (* 4 length-scale)) 'double-float))
	(setf (nth 7 (aref *grid* i j k)) (coerce (- (* j (/ (* 3 length-scale))) (* 4 length-scale)) 'double-float))
	(setf (nth 8 (aref *grid* i j k)) (coerce (- (* k (/ (* 3 length-scale))) (* 4 length-scale)) 'double-float)))))

(defun E-field (ii jj kk)
  "calculate E-field at grid ref ii jj kk from rho at each of the other grid points"
  (declare (type fixnum ii jj kk))
  (let ((ex 0d0) (ey 0d0) (ez 0d0))
    (dotimes (i +grid-points+)
      (dotimes (j +grid-points+)
	(dotimes (k +grid-points+)
 	  (if (nor (< i 0) (>= i (- +grid-points+ 1)) 
		   (< j 0) (>= j (- +grid-points+ 1)) 
		   (< k 0) (>= k (- +grid-points+ 1))
		   (and (= i ii) (= j jj) (= k kk)))
	      (dbind (x y z) (the list (grid-r-vec ii jj kk i j k))
		(let* ((r (q-add x y z))
		       (p/rrr (/ (g-rho i j k) (the double-float r) (the double-float r) (the double-float r))))
		  (incf (the double-float ex) (* (the double-float x) (the double-float p/rrr)))
		  (incf (the double-float ey) (* (the double-float y) (the double-float p/rrr)))
		  (incf (the double-float ez) (* (the double-float z) (the double-float p/rrr)))))))))
    (list ex ey ez)))

(defun B-field (ii jj kk)
  "calculate B-field at grid ref ii jj kk from J at each of the other grid points"
  (declare (type fixnum ii jj kk))
  (let ((bx 0d0) (by 0d0) (bz 0d0))
    (dotimes (i +grid-points+)
      (dotimes (j +grid-points+)
	(dotimes (k +grid-points+)
	  (if (nor (< i 0) (>= i (- +grid-points+ 1)) 
		   (< j 0) (>= j (- +grid-points+ 1)) 
		   (< k 0) (>= k (- +grid-points+ 1))
		   (and (= i ii) (= j jj) (= k kk)))
	      (dbind (x y z) (grid-r-vec ii jj kk i j k)
		(let ((r (q-add x y z)))
		  (incf (the double-float bx) (/ (- (* (g-jy i j k) (the double-float z)) (* (g-jz i j k) (the double-float y))) (the double-float r)(the double-float r)(the double-float r)))
		  (incf (the double-float by) (/ (- (* (g-jz i j k) (the double-float x)) (* (g-jx i j k) (the double-float z))) (the double-float r)(the double-float r)(the double-float r)))
		  (incf (the double-float bz) (/ (- (* (g-jx i j k) (the double-float y)) (* (g-jy i j k) (the double-float x))) (the double-float r)(the double-float r)(the double-float r)))
		  ))))))
    (list bx by bz)))
		  

(defun scalar-pot (x y z &optional (self -1))
  "scalar potential phi, no longer used in calcs but used for electrostatic energy readout"
  (declare (type double-float x y z)
	   (type fixnum self))
  (let ((ans 0d0))
    (dotimes (i +n-ions+)
      (if (= i self)()
	  (decf ans (/ (charge i)
		       (the double-float (q-add (- x (posx i)) (- y (posy i)) (- z (posz i))))))))
    ans))

;(defun vector-pot (ip jp kp &optional (is -1) (js -1) (ks -1))
;  (declare (type double-float x y z)
;	   (type fixnum self))
;  (let ((ax 0d0)(ay 0d0)(az 0d0))
;    	  
;      (if (= i self) ()
;	  (let ((r (the double-float (q-add (- x (posx i)) (- y (posy i)) (- z (posz i))))))
;	    (incf ax (/ (* (charge i) (velx i) ) r))
;	    (incf ay (/ (* (charge i) (vely i) ) r))
;	    (incf az (/ (* (charge i) (velx i) ) r)))))
;	(list ax ay az)))
;
;(defun s-grad (fn)
;  "gradient of 3D scalar field given by FN"
;  #'(lambda (x y z &optional (self -1))
;      (declare (type double-float x y z))
;      (list (funcall (partial-diff fn 0 4) x y z self)
;	    (funcall (partial-diff fn 1 4) x y z self)
;	    (funcall (partial-diff fn 2 4) x y z self))))
;
;(defun s-curl (fn)
;  "curl of 3D vector field given by FN"
;  (let ((d/dx (partial-diff-vec fn 0 4))
;	(d/dy (partial-diff-vec fn 1 4))
;	(d/dz (partial-diff-vec fn 2 4)))
;  #'(lambda (x y z &optional (self -1))
;      (declare (type double-float x y z))
;      (let ((ddx (funcall d/dx x y z self))
;	    (ddy (funcall d/dy x y z self))
;	    (ddz (funcall d/dz x y z self)))
;	(list (- (the double-float (caddr ddy))
;		 (the double-float (cadr ddz)))
;	      (- (the double-float (car ddz))
;		 (the double-float (caddr ddx)))
;	      (- (the double-float (cadr ddx))
;		 (the double-float (car ddy))))))))

(defun q-add (x y z)
  "quad-add optimized for vectors of length 3 with only double-floats"
  (declare (type double-float x y z))
  (sqrt (+ (* x x) (* y y) (* z z))))

(defmemo grid-r (i1 j1 k1 i2 j2 k2)
;the distance between 2 grid points, only +grid-points+³ posible answers so memoize them
  (apply #'q-add (vec- (g-pos i1 j1 k1)
		       (g-pos i2 j2 k2))))

(defmemo grid-r-vec (i1 j1 k1 i2 j2 k2)
;the vector between 2 grid points, only +grid-points+³ posible answers so memoize them
  (vec- (g-pos i2 j2 k2)
	(g-pos i1 j1 k1)))

;(defun update-grid (&optional (ext-B-field #'(lambda (x y z) (declare (ignore x y z)) '(0 0 0))))
;  (dotimes (i +grid-points+)
;    (dotimes (j +grid-points+)
;      (dotimes (k +grid-points+)
;	(let* ((pos (list (coerce (- (* i 2/3) 2) 'double-float) (coerce (- (* j 2/3) 2) 'double-float) (coerce (- (* k 2/3) 2) 'double-float)))
;	       (E-field (s-grad (memoize #'scalar-pot)))
;	       (B-field (s-curl (memoize #'vector-pot)))
;	       (B (vec+ (apply ext-b-field pos) (apply B-field pos)))
;	       (E (apply E-field pos)))
;	  (setf (aref *grid* i j k) (append E B Pos)))))))

(defun update-grid (&optional (ext-B-field #'(lambda (x y z) (declare (ignore x y z)) '(0 0 0))))
  "functon for setting rho and J at each of the grid points then calculating the E and B fields from these new values"
  (dotimes (i +grid-points+) (dotimes (j +grid-points+) (dotimes (k +grid-points+) (setf (g-rho i j k) 0d0
											 (g-jx i j k) 0d0
											 (g-jy i j k) 0d0
											 (g-jz i j k) 0d0))))
  (dotimes (n +n-ions+)
    (dbind (i j k) (list (the double-float (* 1.5d0 (+ 2d0 (posx n)))) 
			 (the double-float(* 1.5d0 (+ 2d0 (posy n)))) 
			 (the double-float(* 1.5d0 (+ 2d0 (posz n))))) ;position of ion in terms of fractional grid places
      (multiple-value-bind (i_a i_b) (floor (the double-float i)) ;split fraction grid place into referance to grid point to the "left" (i.e more negative referance) and fraction between that and next grid
	(multiple-value-bind (j_a j_b) (floor (the double-float j))  ; ie for grid points at 0,0 2,0 0,2 2,2 an ion at 2,1.5 would be i_a=1, i_b=0, j_a=0, j_b=0.75 
	  (multiple-value-bind (k_a k_b) (floor (the double-float k))
	    (if (nor (< i 0) (>= i (- +grid-points+ 1)) (< j 0) (>= j (- +grid-points+ 1)) (< k 0) (>= k (- +grid-points+ 1)))
		(let ((A (+ (/ (q-add i_b j_b k_b))
			    (/ (q-add i_b j_b (- 1d0 k_b)))
			    (/ (q-add i_b (- 1d0 j_b) k_b))
			    (/ (q-add (- 1d0 i_b) j_b k_b))
			    (/ (q-add (- 1d0 i_b) (- 1d0 j_b) k_b))
			    (/ (q-add (- 1d0 i_b) j_b (- 1d0 k_b)))
			    (/ (q-add i_b (- 1d0 j_b) (- 1d0 k_b)))           ; A is the sum of the reciprocals of the distances between the ion and the 8 points at the corners of the cell it is in.
			    (/ (q-add (- 1d0 i_b) (- 1d0 j_b) (- 1d0 k_b))))) ; this ensures when the charge of an ion is distributed between points the total to all of them is rqual to the ions charge
		      (jx (* (velx n) (charge n)))
		      (jy (* (vely n) (charge n)))
		      (jz (* (velz n) (charge n))))
		  
;setting new rho and J using 
; rho_at_point = sum_over_all_ions_in_a_cell_adjacent_to_this_point(charge_of_ion/(A*distance_between_ion_and_point))
;ditto for each componant of J but using current in that direction rather than charge
		  (incf (g-rho      i_a j_a k_a )             (/ (charge n) A (q-add i_b j_b k_b)))
		  (incf (g-rho      i_a j_a (+ 1 k_a) )       (/ (charge n) A (q-add i_b j_b (- 1d0 k_b))))
		  (incf (g-rho      i_a (+ 1 j_a) k_a )       (/ (charge n) A (q-add i_b (- 1d0 j_b) k_b)))
		  (incf (g-rho (+ 1 i_a) j_a k_a )            (/ (charge n) A (q-add (- 1d0 i_b) j_b k_b)))
		  (incf (g-rho (+ 1 i_a) (+ 1 j_a) k_a )      (/ (charge n) A (q-add (- 1d0 i_b) (- 1d0 j_b) k_b)))
		  (incf (g-rho (+ 1 i_a) j_a (+ 1 k_a))       (/ (charge n) A (q-add (- 1d0 i_b) j_b (- 1d0 k_b))))
		  (incf (g-rho      i_a (+ 1 j_a) (+ 1 k_a))  (/ (charge n) A (q-add i_b (- 1d0 j_b) (- 1d0 k_b))))
		  (incf (g-rho (+ 1 i_a) (+ 1 j_a) (+ 1 k_a)) (/ (charge n) A (q-add (- 1d0 i_b) (- 1d0 j_b) (- 1d0 k_b)))) 
		  
		  (incf (g-jx      i_a j_a k_a )             (/ jx A (q-add i_b j_b k_b)))
		  (incf (g-jx      i_a j_a (+ 1 k_a) )       (/ jx A (q-add i_b j_b (- 1d0 k_b))))
		  (incf (g-jx      i_a (+ 1 j_a) k_a )       (/ jx A (q-add i_b (- 1d0 j_b) k_b)))
		  (incf (g-jx (+ 1 i_a) j_a k_a )            (/ jx A (q-add (- 1d0 i_b) j_b k_b)))
		  (incf (g-jx (+ 1 i_a) (+ 1 j_a) k_a )      (/ jx A (q-add (- 1d0 i_b) (- 1d0 j_b) k_b)))
		  (incf (g-jx (+ 1 i_a) j_a (+ 1 k_a))       (/ jx A (q-add (- 1d0 i_b) j_b (- 1d0 k_b))))
		  (incf (g-jx      i_a (+ 1 j_a) (+ 1 k_a))  (/ jx A (q-add i_b (- 1d0 j_b) (- 1d0 k_b))))
		  (incf (g-jx (+ 1 i_a) (+ 1 j_a) (+ 1 k_a)) (/ jx A (q-add (- 1d0 i_b) (- 1d0 j_b) (- 1d0 k_b))))
		  
		  (incf (g-jy      i_a j_a k_a ) (/ jy A (q-add i_b j_b k_b)))
		  (incf (g-jy      i_a j_a (+ 1 k_a) )       (/ jy A (q-add i_b j_b (- 1d0 k_b))))
		  (incf (g-jy      i_a (+ 1 j_a) k_a )       (/ jy A (q-add i_b (- 1d0 j_b) k_b)))
		  (incf (g-jy (+ 1 i_a) j_a k_a )            (/ jy A (q-add (- 1d0 i_b) j_b k_b)))
		  (incf (g-jy (+ 1 i_a) (+ 1 j_a) k_a )      (/ jy A (q-add (- 1d0 i_b) (- 1d0 j_b) k_b)))
		  (incf (g-jy (+ 1 i_a) j_a (+ 1 k_a))       (/ jy A (q-add (- 1d0 i_b) j_b (- 1d0 k_b))))
		  (incf (g-jy      i_a (+ 1 j_a) (+ 1 k_a))  (/ jy A (q-add i_b (- 1d0 j_b) (- 1d0 k_b))))
		  (incf (g-jy (+ 1 i_a) (+ 1 j_a) (+ 1 k_a)) (/ jy A (q-add (- 1d0 i_b) (- 1d0 j_b) (- 1d0 k_b))))
		  
		  (incf (g-jz      i_a j_a k_a )             (/ jz A (q-add i_b j_b k_b)))
		  (incf (g-jz      i_a j_a (+ 1 k_a) )       (/ jz A (q-add i_b j_b (- 1d0 k_b))))
		  (incf (g-jz      i_a (+ 1 j_a) k_a )       (/ jz A (q-add i_b (- 1d0 j_b) k_b)))
		  (incf (g-jz (+ 1 i_a) j_a k_a )            (/ jz A (q-add (- 1d0 i_b) j_b k_b)))
		  (incf (g-jz (+ 1 i_a) (+ 1 j_a) k_a )      (/ jz A (q-add (- 1d0 i_b) (- 1d0 j_b) k_b)))
		  (incf (g-jz (+ 1 i_a) j_a (+ 1 k_a))       (/ jz A (q-add (- 1d0 i_b) j_b (- 1d0 k_b))))
		  (incf (g-jz      i_a (+ 1 j_a) (+ 1 k_a))  (/ jz A (q-add i_b (- 1d0 j_b) (- 1d0 k_b))))
		  (incf (g-jz (+ 1 i_a) (+ 1 j_a) (+ 1 k_a)) (/ jz A (q-add (- 1d0 i_b) (- 1d0 j_b) (- 1d0 k_b)))))))))))
                  
  (dotimes (i +grid-points+)
    (dotimes (j +grid-points+)
      (dotimes (k +grid-points+)
	(dbind ((ex ey ez)(bx by bz)) (list (E-field i j k) (vec+ (B-field i j k) (apply ext-b-field (g-pos i j k))))
	  (setf (g-ex i j k) ex   (g-ey i j k) ey   (g-ez i j k) ez
		(g-bx i j k) bx   (g-by i j k) by   (g-bz i j k) bz)))))) ;use E-field and B-field to calc the new E and B fields at each point and set them

(defun interpolate (pos)
  (dbind (i j k) (mapcar #[the double-float (* 1.5d0 (+ 2d0 (the double-float $)))] pos)
    (multiple-value-bind (i_a i_b) (floor i)    ;i_a is the x reference of the cell in position below the position (i.e the first cell with a position < the position) 
      (multiple-value-bind (j_a j_b) (floor j)  ;i_b is the fraction of the way it it between the inferior grid reference and the superior one
	(multiple-value-bind (k_a k_b) (floor k)
	  (if (or (< i 0) (>= i (- +grid-points+ 1)) (< j 0) (>= j (- +grid-points+ 1)) (< k 0) (>= k (- +grid-points+ 1)))
	      '(0d0 0d0 0d0 0d0 0d0 0d0)
	      (dbind ((ex1 ey1 ez1 bx1 by1 bz1)  (ex2 ey2 ez2 bx2 by2 bz2) ;bind E and B values for the 8 points
		      (ex3 ey3 ez3 bx3 by3 bz3)  (ex4 ey4 ez4 bx4 by4 bz4)
		      (ex5 ey5 ez5 bx5 by5 bz5)  (ex6 ey6 ez6 bx6 by6 bz6)
		      (ex7 ey7 ez7 bx7 by7 bz7)  (ex8 ey8 ez8 bx8 by8 bz8))
		     (list (aref *grid*   i_a       j_a     k_a)     (aref *grid*   i_a       j_a   (1+ k_a))
			   (aref *grid*   i_a    (1+ j_a)   k_a)     (aref *grid* (1+ i_a)    j_a     k_a)
			   (aref *grid* (1+ i_a) (1+ j_a)   k_a)     (aref *grid* (1+ i_a)    j_a   (1+ k_a))
			   (aref *grid*   i_a    (1+ j_a) (1+ k_a))  (aref *grid* (1+ i_a) (1+ j_a) (1+ k_a)))
		(let* ((A (q-add   i_b       j_b     k_b))
		       (b (q-add   i_b       j_b   (1+ k_b)))
		       (c (q-add   i_b    (1+ j_b)   k_b))
		       (d (q-add (1+ i_b)    j_b     k_b))
		       (e (q-add (1+ i_b) (1+ j_b)   k_b))
		       (f (q-add (1+ i_b)    j_b   (1+ k_b)))
		       (g (q-add   i_b    (1+ j_b) (1+ k_b)))
		       (h (q-add (1+ i_b) (1+ j_b) (1+ k_b)))
		       (AA (+ (/ (the double-float a))(/ (the double-float b))(/ (the double-float c))(/ (the double-float d))
			      (/ (the double-float e))(/ (the double-float f))(/ (the double-float g))(/ (the double-float h)))));distance from the ion to each of the 8 points and the sum of their reciprocals (for normilisation)
		  (declare (type double-float a b c d e f g h aa))
		  (list (+ (/ (the double-float ex1) AA A)(/ (the double-float ex2) AA b)(/ (the double-float ex3) AA c)(/ (the double-float ex4) AA d)(/ (the double-float ex5) AA e)(/ (the double-float ex6) AA f)(/ (the double-float ex7) AA g)(/ (the double-float ex8) AA h)) ;E_x
			(+ (/ (the double-float ey1) AA A)(/ (the double-float ey2) AA b)(/ (the double-float ey3) AA c)(/ (the double-float ey4) AA d)(/ (the double-float ey5) AA e)(/ (the double-float ey6) AA f)(/ (the double-float ey7) AA g)(/ (the double-float ey8) AA h)) ;E_y
			(+ (/ (the double-float ez1) AA A)(/ (the double-float ez2) AA b)(/ (the double-float ez3) AA c)(/ (the double-float ez4) AA d)(/ (the double-float ez5) AA e)(/ (the double-float ez6) AA f)(/ (the double-float ez7) AA g)(/ (the double-float ez8) AA h)) ;E_z
			(+ (/ (the double-float bx1) AA A)(/ (the double-float bx2) AA b)(/ (the double-float bx3) AA c)(/ (the double-float bx4) AA d)(/ (the double-float bx5) AA b)(/ (the double-float bx6) AA f)(/ (the double-float bx7) AA g)(/ (the double-float bx8) AA h)) ;B_x
			(+ (/ (the double-float by1) AA A)(/ (the double-float by2) AA b)(/ (the double-float by3) AA c)(/ (the double-float by4) AA d)(/ (the double-float by5) AA b)(/ (the double-float by6) AA f)(/ (the double-float by7) AA g)(/ (the double-float by8) AA h)) ;B_y
			(+ (/ (the double-float bz1) AA A)(/ (the double-float bz2) AA b)(/ (the double-float bz3) AA c)(/ (the double-float bz4) AA d)(/ (the double-float bz5) AA b)(/ (the double-float bz6) AA f)(/ (the double-float bz7) AA g)(/ (the double-float bz8) AA h)) ;B_z
			)))))))))


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




(defun main (&optional (ticks 100) Bf)
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
	  (with-open-file (stream5 "dat/dat5"   :direction :output :if-exists :supersede :if-does-not-exist :create)
	    (with-open-file (stream6 "dat/dat6"   :direction :output :if-exists :supersede :if-does-not-exist :create)
	      (with-open-file (stream7 "dat/dat7"   :direction :output :if-exists :supersede :if-does-not-exist :create)
		(with-open-file (stream8 "dat/dat8"   :direction :output :if-exists :supersede :if-does-not-exist :create)
		  (with-open-file (stream9 "dat/dat9"   :direction :output :if-exists :supersede :if-does-not-exist :create)
		    (with-open-file (streamEE "dat/datEEs" :direction :output :if-exists :supersede :if-does-not-exist :create)
		      (with-open-file (streamEK "dat/datEK"  :direction :output :if-exists :supersede :if-does-not-exist :create)
			
	      
	      (dotimes (i ticks)
		(print i)
		(format streamEK "~&~S " i)(E-kinetic streamEK)
		(format streamEE "~&~S " i)(E-electro streamEE)
		(format stream9 "~&~S ~S ~S~&" (coerce (posx 17) 'single-float) (coerce (posy 17) 'single-float) (coerce (posz 17) 'single-float))
		(format stream8 "~&~S ~S ~S~&" (coerce (posx 15) 'single-float) (coerce (posy 15) 'single-float) (coerce (posz 15) 'single-float))
		(format stream7 "~&~S ~S ~S~&" (coerce (posx 13) 'single-float) (coerce (posy 13) 'single-float) (coerce (posz 13) 'single-float))
		(format stream6 "~&~S ~S ~S~&" (coerce (posx 11) 'single-float) (coerce (posy 11) 'single-float) (coerce (posz 11) 'single-float))
		(format stream5 "~&~S ~S ~S~&" (coerce (posx 9)  'single-float) (coerce (posy 9)  'single-float) (coerce (posz 9)  'single-float))
		(format stream4 "~&~S ~S ~S~&" (coerce (posx 7)  'single-float) (coerce (posy 7)  'single-float) (coerce (posz 7)  'single-float))
		(format stream3 "~&~S ~S ~S~&" (coerce (posx 5)  'single-float) (coerce (posy 3)  'single-float) (coerce (posz 3)  'single-float))
		(format stream2 "~&~S ~S ~S~&" (coerce (posx 3)  'single-float) (coerce (posy 5)  'single-float) (coerce (posz 5)  'single-float))
		(format stream1 "~&~S ~S ~S~&" (coerce (posx 1)  'single-float) (coerce (posy 1)  'single-float) (coerce (posz 1)  'single-float))
    
    
		(if bf (update-grid bf) (update-grid))
		(dotimes (n +n-ions+)
		  (incf (posx n) (* (velx n) delta-t))
		  (incf (posy n) (* (vely n) delta-t))
		  (incf (posz n) (* (velz n) delta-t))
		  (setf (posz n) (- (mod (+ 2d0 (posz n)) 4d0) 2d0))
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
		      (incf (velz n) (* 0.5d0 (the double-float ez) delta-t (/ (charge n) (mass n)))))))))))))))))))

  (with-open-file (stream "dat/dat-final-pos" :direction :output :if-exists :supersede :if-does-not-exist :create)
    (positions stream))
  (with-open-file (stream "dat/dat-final-vel" :direction :output :if-exists :supersede :if-does-not-exist :create)
    (velocities stream))
  (with-open-file (stream "dat/dat-final-x-vx" :direction :output :if-exists :supersede :if-does-not-exist :create)
    (x-v-vx stream)))

