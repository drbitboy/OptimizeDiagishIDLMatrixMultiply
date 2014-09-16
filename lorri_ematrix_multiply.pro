;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Mutiply an NxN matrix by the epsilon matrix (EMatrix)
;;;
;;; - N.B. N must be less than or equal to the size of EMatrix
;;;
;;; - Provides about a 10x speedup over NxN # EMatrix for N=1024
;;;   - LORRI desmear calibration code does four such multiplies
;;;     - 1 billion multiplies and adds each
;;;     - each takes about 0.5s for N=1024
;;;
;;; - For a random 10241024 input matrix (=randomu(Seed,1024,1024,/dou),
;;;   the fractional errors are less than 1E-14
;;;
;;; - Required form of EMatrix (from SOC_INST_ICD):
;;;
;;;                   /  Tf1/Tavg if k>j
;;;   EMatrix[k,j] = <          1 if k=j
;;;                   \  Tf2/Tavg if k<j
;;;
;;;   Tavg = (Tf1 + Tf2) / 2
;;;
;;; - As of 2014-07-19:  Tf1/Tavg = 1.044; Tf2/Tavg = 0.956
;;;
;;; 2014-07-19 BTCarcich Initial version
;;; 2014-09-16 BTCarcich Better version, ~10x speedup wrt #
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Details:
;;;
;;;   [K] = [IMG] # [E]
;;;
;;; Since that multiplies columns of [IMG] by rows of [E], there is no
;;; interaction between columns of [IMG], so for analysis we only need
;;; to look at one column of [IMG], the dot (inner) products of which,
;;; with the rows of [E], become one column of the result [K]. Expanding
;;; that idea to one column of a 3x3 matrix with elements x, y and z:
;;;
;;;   [Kx]   [x]   [1    1+d  1+d]
;;;   [Ky] = [y] # [1-d  1    1+d]
;;;   [Kz] = [z]   [1-d  1-d  1  ]
;;;
;;; Looking specifically at element Ky above (one element of [K],
;;; corresponding to y in [IMG], broken down to a formula using only
;;; scalars):
;;;
;;;   Ky = x * (1-d)  +  y * 1  + z * (1+d)
;;;
;;; Generalizing for any y in a column of any length:
;;;
;;;   Ky = sumx * (1-d)  +  y * 1  + sumz * (1+d)
;;;
;;; Where scalars sumx and sumz and sums of all values above and below,
;;; respectively, y in that column of [IMG]. N.B. sumx and sumz are
;;; values specific to element y.
;;;
;;; Rearranging:
;;;
;;;   Ky = (sumx + y + sumz) - d * (sumx - sumz)
;;;
;;;   Ky = tot - d * (sumx - sumz)
;;;
;;; where
;;;
;;;   tot = (sumx + y + sumz)
;;;
;;; I.e. tot is the sum of all values in the column
;;;      (e.g. in IDL: tot = total(IMG,2)).
;;;
;;; So to this point I have basically duplicated Diego's work; the rest
;;; of this analysis converts that last equation for Ky into a form
;;; suitable for speedy evaluation in IDL.
;;;
;;; Solving the tot equation for sumz:
;;;
;;;   sumz = tot - (y + sumx)
;;;
;;; Substituting back into Ky:
;;;
;;;   Ky = tot - (sumx - (tot - (y + sumx)))
;;;
;;;   Ky = tot - ((2 * sumx) + y - tot)
;;;
;;;   Ky = tot + (tot - ((2 * sumx) + y)
;;;
;;; Using sumxy to represent the sum of all values in the column from
;;; the top down to, and including, y
;;; (IDL: [SUMXY] = total([IMG],2,/CUMULATIVE))
;;;
;;;   sumxy = sumx + y
;;;
;;; and
;;;
;;;   sumx = sumxy - y
;;;
;;; Substituting back into Ky:
;;;
;;;   Ky = tot + (tot - ((2 * (sumxy - y)) + y)
;;;
;;;   Ky = tot + (tot + y - (2 * sumxy))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION LORRI_EMATRIX_MULTIPLY, NxN, EMatrix, delta=deltaArg
  ;;; Get delta from delta arugment if present, else from EMatrix[1,0]
  DELTA = N_ELEMENTS(deltaArg) eq 1L $
          ? DOUBLE(deltaArg[0])      $
          : (EMatrix[1,0] - 1)
  ;;; Get number of rows in array
  NROWS = (SIZE(NxN,/DIMENSIONS))[1]
  ;;; Calculate arra of cumulative sums down columns
  SUMXY = TOTAL(NxN,2,/CUMULATIVE)
  ;;; Duplicate column totals (last row of SUMXY) into complete matrix
  TOT   = SUMXY[*,NROWS-1] # REPLICATE(1,NROWS,1d0)
  ;;; Return array equivalent to (NxN # EMatrix)
  return, TOT + (DELTA * (TOT + NxN - (2 * SUMXY)))
  ;;; N.B. Multiplying by integer 2 appears to be faster than by 2D0
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Test code
;;;
;;; Usage:
;;;
;;;  echo .r lorri_ematrix_multiply.pro | TEST_LORRI_EMATRIX_MULTIPLY=yes idl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if getenv("TEST_LORRI_EMATRIX_MULTIPLY") ne '' then begin
  !quiet=1b

  if n_elements(NxN) eq 0L then NxN = randomu(iseed,1024,1024,/double)
  EMatrix = double(readfits(getenv('HOME') + '/pipeline/level2/lor/cal/dsmear_ematrix_1024.fit'))
  dims = size(NxN,/dim)
  nrows = dims[1]
  EMatrix = EMatrix[0:nrows-1,0:nrows-1]

  if n_elements(nTest) ne 1L then nTest = 10L

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  print,'Starting one'
  t0=systime(1,/seconds)
  for i=0L,nTest-1L do multone = NxN#EMatrix
  print,f='(f10.2,3x,a)',systime(1,/seconds)-t0,'Finished one'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  print,'Starting lorri_ematrix_multiply'
  t0=systime(1,/seconds)
  for i=0L,nTest-1L do multlorri_ematrix_multiply = lorri_ematrix_multiply(NxN,EMatrix)
  print,f='(f10.2,3x,a)',systime(1,/seconds)-t0,'Finished lorri_ematrix_multiply'

  help,max(abs(multone-multlorri_ematrix_multiply))
  help,max(abs(multone-multlorri_ematrix_multiply)/abs(multone))

  help,where(abs(multone-multlorri_ematrix_multiply)/abs(multone) gt 1d-15)
  help,where(abs(multone-multlorri_ematrix_multiply)/abs(multone) gt 1d-14)
  help,where(abs(multone-multlorri_ematrix_multiply)/abs(multone) gt 1d-13)
  help,where(abs(multone-multlorri_ematrix_multiply)/abs(multone) gt 1d-12)

endif

end
