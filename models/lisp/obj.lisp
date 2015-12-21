;;; This is a simple script to strip an obj to the bare minimum (vertices + faces), in
;;; order to be ready for processing by a crude obj file reader. It is built around the
;;; notion of a line processor, which receives an input line as an argument and returns a
;;; list of output lines as its value.

(ql:quickload :split-sequence)

(defun strip (input-pathname)
  (process-obj input-pathname
	       (make-output-pathname input-pathname
				     'strip)
	       #'strip-line))

(defun scale (input-pathname &key (x 1) (y 1) (z 1))
  (process-obj input-pathname
	       (make-output-pathname input-pathname
				     'scale)
	       (make-scaling-line-processor x y z)))

(defun make-output-pathname (input-pathname operation)
  (make-pathname :name (concatenate 'string
				    (pathname-name input-pathname)
				    "_"
				    (string-downcase (symbol-name operation)))
		 :defaults input-pathname))

(defun process-obj (input-pathname output-pathname line-processor)
  (with-open-file (in input-pathname
		      :direction :input)
    (with-open-file (out output-pathname
			 :direction :output
			 :if-does-not-exist :create
			 :if-exists :supersede)
      (loop (let ((line (read-line in nil nil)))
	      (if line
		  (format out "~{~a~%~}" (funcall line-processor line))
		  (return)))))))

(defun strip-line (line)
  (let ((tokens (split-sequence:split-sequence #\space
					       line
					       :remove-empty-subseqs t)))
    (cond ((string= (car tokens) "v")
	   (list line))
	  ((string= (car tokens) "f")
	   (strip-face-tokens (cdr tokens)))
	  (t nil))))

(defun strip-face-tokens (tokens)
  (let* ((vertex-indices (mapcar (lambda (token)
				   (car (split-sequence:split-sequence #\/ token)))
				 tokens))
	 (triangle-count (- (length vertex-indices) 2))
	 (result-lines nil))
    (dotimes (i triangle-count result-lines)
      (push (format nil
		    "f ~a ~a ~a"
		    (nth 0 vertex-indices)
		    (nth (+ 1 i) vertex-indices)
		    (nth (+ 2 i) vertex-indices))
	    result-lines))))

(defun make-scaling-line-processor (x y z)
  (let ((factors (list x y z)))
    (lambda (line)
      (let ((tokens (split-sequence:split-sequence #\space
						   line
						   :remove-empty-subseqs t)))
	(cond ((string= (car tokens) "v")
	       (list (apply #'format
			    nil
			    "v ~a ~a ~a"
			    (mapcar (lambda (input-token factor)
				      (* factor
					 (with-input-from-string (in input-token)
					   (read in nil nil))))
				    (cdr tokens)
				    factors))))
	      (t (list line)))))))
