;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Mutiply an NxN matrix by the epsilon matrix (EMatrix)
;;;
;;; - N.B. N must be less than or equal to the size of EMatrix
;;;
;;; - Provides about a 6x speedup over NxN # EMatrix for N=1024
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
;;; 2014-09-19 BTCarcich
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function lorri_ematrix_multiply,NxN,EMatrix,NxNplus=NxNplus,NxNminus=NxNminus
common two_two_cmn, onerow1k, onerowQk $
                  , iw1kRotate7, iwQkRotate7 $
                  , iw1kShift01, iwQkShift01 $
                  , iw1kRotate7Shift01, iwQkRotate7Shift01 $
                  , commonIsBuilt

    ;;; Save indices to perform rotation and shifting

    if ~ keyword_set(commonIsBuilt) then begin
      iw1k = lindgen(1024,1024)
      iwQk = lindgen(256,256)

      onerow1k=iw1k[0:1023]
      onerowQk=iwQk[0:255]

      iw1kRotate7 = rotate(iw1k,7)
      iwQkRotate7 = rotate(iwQk,7)

      iw1kShift01 = shift(iw1k,0,1)
      iwQkShift01 = shift(iwQk,0,1)

      iw1kRotate7Shift01 = iw1kRotate7[iw1kShift01]
      iwQkRotate7Shift01 = iwQkRotate7[iwQkShift01]

      commonIsBuilt = 1b
    endif
    nrows = (size(NxN,/dim))[1]

    ;;; Get indices specific to NxN

    case nrows of
    1024: begin
            onerow = onerow1k
            iwRotate7 = iw1kRotate7
            iwShift01 = iw1kShift01
            iwRotate7Shift01 = iw1kRotate7Shift01
          end
     256: begin
            onerow = onerowQk
            iwRotate7 = iwQkRotate7
            iwShift01 = iwQkShift01
            iwRotate7Shift01 = iwQkRotate7Shift01
          end
    else: begin
            iwNRows = lindgen(nRows,nRows)
            onerow = iwNRows[0:nRows-1]
            iwRotate7 = rotate(iwNRows,7)
            iwShift01 = shift(iwNRows,0,1)
            iwRotate7Shift01 = iwRotate7[iwShift01]
          end
    endcase

    ;;; Expand NxN input into two additional matrices:
    ;;;
    ;;; - NxNplus will be mutiplied by Tf1/Tavg
    ;;;
    ;;; - NxNminus will be mutiplied by Tf2/Tavg

    NxNplus = NxN[iwRotate7Shift01]
    NxNplus[onerow] = 0d0
    NxNplus = total(/double,/cumul,temporary(NxNplus),2)
    NxNplus = NxNplus[iwRotate7]

    NxNminus = NxN[iwShift01]
    NxNminus[onerow] = 0d0
    NxNminus = total(/double,/cumul,temporary(NxNminus),2)

    ;;; Perform the equivalent of 
    return, NxNplus*EMatrix[1,0] + NxN + NxNminus*EMatrix[0,1]

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Backup if READFITS.PRO is not in path
function readfits_backup
  print,'### USING READFITS_BACKUP'
  delta = 0.044d0
  sz = 1024L
  rtn = make_array(sz,sz,value=1d0)
  iw1Kx1K = lindgen(sz,sz)
  cols = iw1Kx1K MOD sz
  rows = iw1Kx1K / sz
  rtn[where( cols gt rows)] += delta
  rtn[where( cols lt rows)] -= delta
  return,rtn
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Test code
;;;
;;; Usage:
;;;
;;;  echo .r lorri_ematrix_multiply.pro | TEST_LORRI_EMATRIX_MULTIPLY=yes idl
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if getenv("TEST_LORRI_EMATRIX_MULTIPLY") ne '' $
or n_elements(nTest) eq 1L then begin
  !quiet=1b

  if n_elements(NxN) eq 0L then NxN = randomu(iseed,1024,1024,/double)

  catcherr = 0L
  catch,catcherr
  if catcherr eq 0L then begin
    EMatrix = double(readfits('dsmear_ematrix_1024.fit'))
  endif else begin
    EMatrix = readfits_backup()
  endelse
  catch,/cancel

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
  for i=0L,nTest-1L do multlorri_ematrix_multiply = lorri_ematrix_multiply(NxN,EMatrix,NxNplus=NxNplus,NxNminus=NxNminus)
  print,f='(f10.2,3x,a)',systime(1,/seconds)-t0,'Finished lorri_ematrix_multiply'

  help,where(abs(multone-multlorri_ematrix_multiply)/abs(multone) gt 1d-15)
  help,where(abs(multone-multlorri_ematrix_multiply)/abs(multone) gt 1d-14)
  help,where(abs(multone-multlorri_ematrix_multiply)/abs(multone) gt 1d-13)
  help,where(abs(multone-multlorri_ematrix_multiply)/abs(multone) gt 1d-12)

endif

end
