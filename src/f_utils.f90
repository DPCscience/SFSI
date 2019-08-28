! ----------------------------------------------------------
! Save a binary file in a Fortran format
! ----------------------------------------------------------

SUBROUTINE writeBinFile2(filename,nfilename,nrows,ncols,sizevar,X,ios)
IMPLICIT NONE
INTEGER*4, intent(in) :: nrows, ncols, sizevar, nfilename
INTEGER, intent(out) :: ios
INTEGER*4 :: filename(nfilename)
CHARACTER*100 :: filename2
REAL*8, TARGET :: X(nrows,ncols)
REAL*8, POINTER :: pX(:,:)
INTEGER :: i

pX => X

filename2=ACHAR(filename(1))
DO i = 1,(nfilename-1)
filename2=filename2(1:i)//ACHAR(filename(i+1))
END DO

OPEN(UNIT=11, FILE=filename2, STATUS="NEW", ACCESS="STREAM",IOSTAT=ios)
WRITE(11) nrows
WRITE(11) ncols
WRITE(11) sizevar

DO i = 1,nrows
IF (sizevar .EQ. 4) THEN
WRITE(11) real(pX(i,:),4)
ELSE
WRITE(11) pX(i,:)
END IF
END DO

CLOSE(UNIT=11)

END SUBROUTINE


! ----------------------------------------------------------
! Read a binary file in a Fortran format
! ----------------------------------------------------------
SUBROUTINE readBinFile2(filename,nfilename,nsetRow,nsetCol,setRow,setCol,ios,ncols,sizevar,n,p,X)
IMPLICIT NONE
INTEGER*4, intent(in) :: nfilename, nsetRow, nsetCol, ncols, sizevar, n, p
INTEGER*4 :: filename(nfilename), setRow(nsetRow), setCol(nsetCol)
INTEGER*4 :: i, ios
CHARACTER*100 :: filename2
REAL*8, intent(out), TARGET :: X(n,p)
REAL*8, POINTER :: pX(:,:)
REAL*4 lineSingle(ncols)
REAL*8 lineDouble(ncols)

pX => X

filename2=ACHAR(filename(1))
DO i = 1,(nfilename-1)
filename2=filename2(1:i)//ACHAR(filename(i+1))
END DO

OPEN(UNIT=42, FILE=filename2, STATUS="OLD", ACCESS="STREAM",form='unformatted',IOSTAT=ios)
! Skip top 3 lines
CALL FSEEK(42, 12, 0)

! Read lines
DO i = 1,n
IF (nsetRow .GT. 0) THEN
CALL FSEEK(42, 12+ncols*sizevar*(setRow(i)-1), 0)
END IF

IF (sizevar .EQ. 4) THEN
READ(42) lineSingle

IF (nsetCol .GT. 0) THEN
pX(i,:)=lineSingle(setCol)
ELSE
pX(i,:)=lineSingle
END IF

ELSE
READ(42) lineDouble

IF (nsetCol .GT. 0) THEN
pX(i,:)=lineDouble(setCol)
ELSE
pX(i,:)=lineDouble
END IF
END IF
END DO

CLOSE(UNIT=42)

END SUBROUTINE
