(ql:quickload '(:split-sequence))
(use-package :split-sequence)

(defstruct (impulse (:constructor make-impulse (magnitude x y z)))
  (magnitude 0.0 :type short-float)
  (x 0.0 :type short-float)
  (y 0.0 :type short-float)
  (z 0.0 :type short-float))

(defun parse-single-float (string)
  (with-input-from-string (in string)
    (coerce (read in) 'single-float)))

(defun read-impulses (pathname-regex)
  (let* ((files (directory pathname-regex))
	 (frame-count (length files))
	 (impulses (make-list frame-count :initial-element nil))
	 (current-frame 0))
    (dolist (file files impulses)
      (setf (nth current-frame impulses)
	    (with-open-file (in file :direction :input)
	      (read-line in)
	      (let ((frame-impulses nil))
		(loop (let ((impulse-tokens (split-sequence #\,
							    (read-line in nil)
							    :remove-empty-subseqs t)))
			(if impulse-tokens
			    (push (apply #'make-impulse
					 (mapcar #'parse-single-float
						 impulse-tokens))
				  frame-impulses)
			    (return))))
		frame-impulses)))
      (print (incf current-frame)))))
