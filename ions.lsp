(load "/home/tom/programs/lisp/maths.fasl")
(proclaim '(inline scalar-pot vector-pot q-add))
(proclaim '(optimize (speed 3) (compilation-speed 0) (debug 3)))

(defconstant epsi0 8.85418782d-12)
(defconstant mu0  (/ 1 c0 c0 epsi0))
(defconstant e-charge 1.60217646d-19)
(defconstant p-mass 1.67262158d-27)
(defconstant e-mass 9.10938188d-31)

(defparameter *arr-length* 1000)
(defvar *ion-array*
  (make-array `(,*arr-length* 8) :element-type 'double-float))
(defvar *temp* 
  (make-array `(,*arr-length* 8) :element-type 'double-float))

(defmacro posx (n)   `(the double-float (aref *ion-array* ,n 0)))
(defmacro posy (n)   `(the double-float (aref *ion-array* ,n 1)))
(defmacro posz (n)   `(the double-float (aref *ion-array* ,n 2)))
(defmacro momx (n)   `(the double-float (aref *ion-array* ,n 3)))
(defmacro momy (n)   `(the double-float (aref *ion-array* ,n 4)))
(defmacro momz (n)   `(the double-float (aref *ion-array* ,n 5)))
(defmacro charge (n) `(the double-float (aref *ion-array* ,n 6)))
(defmacro mass (n)   `(the double-float (aref *ion-array* ,n 7)))

(defmacro tposx (n)   `(the double-float (aref *temp* ,n 0)))
(defmacro tposy (n)   `(the double-float (aref *temp* ,n 1)))
(defmacro tposz (n)   `(the double-float (aref *temp* ,n 2)))
(defmacro tmomx (n)   `(the double-float (aref *temp* ,n 3)))
(defmacro tmomy (n)   `(the double-float (aref *temp* ,n 4)))
(defmacro tmomz (n)   `(the double-float (aref *temp* ,n 5)))
(defmacro tcharge (n) `(the double-float (aref *temp* ,n 6)))
(defmacro tmass (n)   `(the double-float (aref *temp* ,n 7)))

(dotimes (i *arr-length*)
  (let ((r (gauss-random 2d0 5d0))
	(theta (random (* 2 pi)))
	(z (gauss-random 1d0))
	(charge (* e-charge (expt -1d0 (random 2)))))
    (setf (posx i) (* r (cos theta))
	  (posy i) (* r (sin theta))
	  (posz i) z
	  (momx i) 0d0 (momy i) 0d0 (momz i) 0d0
	  (charge i) charge
	  (mass i) (if (= charge e-charge) e-mass p-mass))))

(defun update-arr ()
  (dotimes (i *arr-length*)
    (setf (posx i) (+ (posx i) (/ (* (tmomx i) 1d-3) (mass i)));each step is 1ms
	  (posy i) (+ (posy i) (/ (* (tmomy i) 1d-3) (mass i)))
	  (posz i) (+ (posz i) (/ (* (tmomz i) 1d-3) (mass i)))
	  (momx i) (tmomx i)
	  (momy i) (tmomy i)
	  (momz i) (tmomz i))
    (setf (tposx i) (posx i)
	  (tposy i) (posy i)
	  (tposz i) (posz i))))

(defun q-add (x y z)
  (declare (type double-float x y z))
  (sqrt (+ (* x x) (* y y) (* z z))))

(defun scalar-pot (x y z &optional (self -1))
  (declare (type double-float x y z)
	   (type fixnum self))
  (let ((ans 0d0))
    (dotimes (i *arr-length*)
      (if (= i self)()
	  (decf ans (/ (charge i)
		       (the double-float (q-add (- x (posx i)) (- y (posy i)) (- z (posz i))))
		       (* 4 pi epsi0)))))
    	ans))

(defun vector-pot (x y z &optional (self -1))
  (declare (type double-float x y z)
	   (type fixnum self))
  (let ((ax 0d0)(ay 0d0)(az 0d0))
    (dotimes (i *arr-length*)
      (if (= i self) ()
	  (let ((r (the double-float (q-add (- x (posx i)) (- y (posy i)) (- z (posz i))))))
	    (incf ax (/ (* (charge i) (momx i) (/ mu0 4 pi)) (mass i) r))
	    (incf ay (/ (* (charge i) (momy i) (/ mu0 4 pi)) (mass i) r))
	    (incf az (/ (* (charge i) (momx i) (/ mu0 4 pi)) (mass i) r)))))
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

(defun toriodial-field (mag r-max width)
  "creates a closure which when given x y z co-ordintates gives the vector of a defined torrodal field

mag is the strength of the toroidal field at its maximum
r-max is the radius in the xy plane at which the field is at maximum
width is the distance from the maximum at which the field has dropped to 1/e of its max value"
  (declare (type double-float mag r-max width))
  #'(lambda (x y z)
      (declare (type double-float x y z))
      (let ((size (* mag  
		     (expt e (- (/ (+ (* z z) 
				      (expt (- (sqrt (+ (* x x) (* y y))) r-max) 2))
				    width width))))))
	(declare (type double-float size)) 
	(cond ((and (= y 0) (> x 0)) (list 0 size 0))    ;
	      ((and (= x 0) (> y 0)) (list (- size) 0 0));
	      ((and (= y 0) (< x 0)) (list 0 (- size) 0));cases where atan theta = n*pi/2
	      ((and (= x 0) (< y 0)) (list size 0 0))    ;
	      ((and (= x 0) (= y 0)) '( 0 0 0))          ;
	      (t (let ((theta (cond ((and (> x 0) (> y 0)) (atan (/ y x)))       ;
				    ((and (< x 0) (> y 0)) (+ pi (atan (/ y x))));general case for theta
				    ((and (< x 0) (< y 0)) (+ pi (atan (/ y x))));
				    ((and (> x 0) (< y 0)) (atan (/ y x))))))    ;
		   ;(declare (type double-float theta))
		   (list (* size (- (sin theta))) (* size (cos theta)) 0)))))))

(defun tick (&optional (Ext-B-Field-fn #'(lambda (x y z)(declare (ignore x y z)) '(0 0 0))))
  "main updating loop calculating the fields at each particles position then storing the new momenta in a separate array then updating the main array"
  (let ((E-f (s-grad (memoize #'scalar-pot)))  ;scalar-pot and vector-pot are not truly functional as they depend upon the state of the main array which isnt a parameter
	(B-f (s-curl (memoize #'vector-pot)))) ;however between updates the array doesnt change so these E-f and B-f can be memoized as they are purely functional
    (dotimes (i *arr-length*)
      (let ((Ef (funcall E-f (posx i) (posy i) (posz i) i))
	    (Bf (vec+ (vec-cross (list (/ (momx i) (mass i))
				       (/ (momy i) (mass i))
				       (/ (momz i) (mass i)))
				 (funcall B-f (posx i) (posy i) (posz i) i))
		      (funcall Ext-B-field-fn (posx i) (posy i) (posz i)))))
	(incf (tmomx i) (* (charge i) 1d-3 (+ (the double-float (car   Ef)) (the double-float (car   Bf)))))
	(incf (tmomy i) (* (charge i) 1d-3 (+ (the double-float (cadr  Ef)) (the double-float (cadr  Bf)))))
	(incf (tmomz i) (* (charge i) 1d-3 (+ (the double-float (caddr Ef)) (the double-float (caddr Bf))))))))
  (update-arr))
