(ql:quickload :split-sequence)
(use-package :split-sequence)

(defun strip (input-pathname)
  "Strip the INPUT-PATHNAME obj file to just vertices and triangular faces."
  (process-obj input-pathname
	       (make-output-pathname input-pathname
				     'strip)
	       #'strip-line))

(defun scale (input-pathname &key (x 1) (y 1) (z 1))
  "Scale INPUT-PATHNAME obj file vertices by factors of X, Y, Z in each dimension
respectively."
  (process-obj input-pathname
	       (make-output-pathname input-pathname
				     'scale)
	       (make-scaling-line-processor x y z)))

(defun make-output-pathname (input-pathname operation)
  "Construct the output pathname by appending the OPERATION name to the filename of the
INPUT-PATHNAME."
  (make-pathname :name (concatenate 'string
				    (pathname-name input-pathname)
				    "_"
				    (string-downcase (symbol-name operation)))
		 :defaults input-pathname))

(defun process-obj (input-pathname output-pathname line-processor)
  "Core helper function, does the file processing from INPUT-PATHNAME to OUTPUT-PATHNAME,
using the LINE-PROCESSOR function supplied, which returns a list of processed output
lines."
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
  "Line processor function. Returns LINE intact if it's a vertex line, performs splitting
into triangular faces for a face LINE, else discards it."
  (let ((tokens (split-sequence #\space
				line
				:remove-empty-subseqs t)))
    (cond ((string= (car tokens) "v")
	   (list line))
	  ((string= (car tokens) "f")
	   (strip-face-tokens (cdr tokens)))
	  (t nil))))

(defun strip-face-tokens (tokens)
  "STRIP-LINE helper creating triangular face lines from one multilateral face line's
tokens."
  (let* ((vertex-indices (mapcar (lambda (token)
				   (car (split-sequence #\/ token)))
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
  "Function which returns a scaling line processor which scales a vertex line's
coordinates by factors of X, Y and Z."
  (let ((factors (list x y z)))
    (lambda (line)
      (let ((tokens (split-sequence #\space
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
