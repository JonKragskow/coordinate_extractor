PROGRAM coordinate_extractor
IMPLICIT NONE

CHARACTER(LEN = 500) :: InputFile, Line, CDummy, OutputFile, FileType
CHARACTER(LEN = 5), ALLOCATABLE :: AtomName(:)
INTEGER :: NAtoms, Idummy, J, AllocationFlag, AC, SOCounter
INTEGER, ALLOCATABLE :: AtomicNumber(:)
REAL(KIND = 8), ALLOCATABLE :: XYZ(:,:)


CALL GET_COMMAND_ARGUMENT(1,InputFile)

IF (InputFile == '-h' .OR. InputFile == '') THEN
    WRITE(6,*) 'Reads coordinates from a Gaussian log or molcas .out/.in file'
    WRITE(6,*)
    WRITE(6,*) 'coordinate_extractor file'
    STOP
END IF

IF (Cut_Past_Char(InputFile,'.') == 'log') THEN
    FileType = 'Gaussian'
ELSEIF (Cut_Past_Char(InputFile,'.') == 'out' .OR. Cut_Past_Char(InputFile,'.') == 'output' .OR. Cut_Past_Char(InputFile,'.') == 'in' .OR. Cut_Past_Char(InputFile,'.') == 'Out' .OR. Cut_Past_Char(InputFile,'.') == 'Output') THEN
    FileType = 'Molcas'
ELSE
    WRITE(6,*) 'Unknown input file type'
    STOP
END IF


OPEN(33, FILE = TRIM(ADJUSTL(InputFile)), STATUS = 'OLD')

IF (FileType == 'Gaussian') THEN

    AllocationFlag = 0
        DO WHILE (.true.)
            READ(33,'(A)', END=98) Line
                Line = ADJUSTL(Line)
                IF (Line(1:7) == 'NAtoms=') THEN
                    READ(Line,*) CDummy, NAtoms
                    IF (AllocationFlag == 0) THEN
                        ALLOCATE(AtomicNumber(NAtoms), XYZ(NAtoms,3), AtomName(NAtoms))
                        AtomicNumber = 0
                        XYZ = 0.0_8
                        AllocationFlag = 1
                    END IF
                END IF
        
                IF (Line(1:21) == 'Standard orientation:') THEN
                    READ(33,*)
                    READ(33,*)
                    READ(33,*)
                    READ(33,*)
                    DO J = 1,NAtoms
                        READ(33,*) Idummy, AtomicNumber(J), Idummy, XYZ(J,1), XYZ(J,2), XYZ(J,3)
                    END DO
                    SOCounter = SOCounter + 1
                END IF
        END DO 
        98 CONTINUE
        
        CLOSE(33)
        
        !Replace Atomic Numbers with Letters in new array
        
        DO J = 1, NAtoms
        
        IF (AtomicNumber(J) == 1) AtomName(J) = 'H'
        IF (AtomicNumber(J) == 2) AtomName(J) = 'He'
        IF (AtomicNumber(J) == 3) AtomName(J) = 'Li'
        IF (AtomicNumber(J) == 4) AtomName(J) = 'Be'
        IF (AtomicNumber(J) == 5) AtomName(J) = 'B'
        IF (AtomicNumber(J) == 6) AtomName(J) = 'C'
        IF (AtomicNumber(J) == 7) AtomName(J) = 'N'
        IF (AtomicNumber(J) == 8) AtomName(J) = 'O'
        IF (AtomicNumber(J) == 9) AtomName(J) = 'F'
        IF (AtomicNumber(J) == 10) AtomName(J) = 'Ne'
        IF (AtomicNumber(J) == 11) AtomName(J) = 'Na'
        IF (AtomicNumber(J) == 12) AtomName(J) = 'Mg'
        IF (AtomicNumber(J) == 13) AtomName(J) = 'Al'
        IF (AtomicNumber(J) == 14) AtomName(J) = 'Si'
        IF (AtomicNumber(J) == 15) AtomName(J) = 'P'
        IF (AtomicNumber(J) == 16) AtomName(J) = 'S'
        IF (AtomicNumber(J) == 17) AtomName(J) = 'Cl'
        IF (AtomicNumber(J) == 18) AtomName(J) = 'Ar'
        IF (AtomicNumber(J) == 19) AtomName(J) = 'K'
        IF (AtomicNumber(J) == 20) AtomName(J) = 'Ca'
        IF (AtomicNumber(J) == 21) AtomName(J) = 'Sc'
        IF (AtomicNumber(J) == 22) AtomName(J) = 'Ti'
        IF (AtomicNumber(J) == 23) AtomName(J) = 'V'
        IF (AtomicNumber(J) == 24) AtomName(J) = 'Cr'
        IF (AtomicNumber(J) == 25) AtomName(J) = 'Mn'
        IF (AtomicNumber(J) == 26) AtomName(J) = 'Fe'
        IF (AtomicNumber(J) == 27) AtomName(J) = 'Co'
        IF (AtomicNumber(J) == 28) AtomName(J) = 'Ni'
        IF (AtomicNumber(J) == 29) AtomName(J) = 'Cu'
        IF (AtomicNumber(J) == 30) AtomName(J) = 'Zn'
        IF (AtomicNumber(J) == 31) AtomName(J) = 'Ga'
        IF (AtomicNumber(J) == 32) AtomName(J) = 'Ge'
        IF (AtomicNumber(J) == 33) AtomName(J) = 'As'
        IF (AtomicNumber(J) == 34) AtomName(J) = 'Se'
        IF (AtomicNumber(J) == 35) AtomName(J) = 'Br'
        IF (AtomicNumber(J) == 36) AtomName(J) = 'Kr'
        IF (AtomicNumber(J) == 37) AtomName(J) = 'Rb'
        IF (AtomicNumber(J) == 38) AtomName(J) = 'Sr'
        IF (AtomicNumber(J) == 39) AtomName(J) = 'Y'
        IF (AtomicNumber(J) == 40) AtomName(J) = 'Zr'
        IF (AtomicNumber(J) == 41) AtomName(J) = 'Nb'
        IF (AtomicNumber(J) == 42) AtomName(J) = 'Mo'
        IF (AtomicNumber(J) == 43) AtomName(J) = 'Tc'
        IF (AtomicNumber(J) == 44) AtomName(J) = 'Ru'
        IF (AtomicNumber(J) == 45) AtomName(J) = 'Rh'
        IF (AtomicNumber(J) == 46) AtomName(J) = 'Pd'
        IF (AtomicNumber(J) == 47) AtomName(J) = 'Ag'
        IF (AtomicNumber(J) == 48) AtomName(J) = 'Cd'
        IF (AtomicNumber(J) == 49) AtomName(J) = 'In'
        IF (AtomicNumber(J) == 50) AtomName(J) = 'Sn'
        IF (AtomicNumber(J) == 51) AtomName(J) = 'Sb'
        IF (AtomicNumber(J) == 52) AtomName(J) = 'Te'
        IF (AtomicNumber(J) == 53) AtomName(J) = 'I'
        IF (AtomicNumber(J) == 54) AtomName(J) = 'Xe'
        IF (AtomicNumber(J) == 55) AtomName(J) = 'Cs'
        IF (AtomicNumber(J) == 56) AtomName(J) = 'Ba'
        IF (AtomicNumber(J) == 57) AtomName(J) = 'La'
        IF (AtomicNumber(J) == 58) AtomName(J) = 'Ce'
        IF (AtomicNumber(J) == 59) AtomName(J) = 'Pr'
        IF (AtomicNumber(J) == 60) AtomName(J) = 'Nd'
        IF (AtomicNumber(J) == 61) AtomName(J) = 'Pm'
        IF (AtomicNumber(J) == 62) AtomName(J) = 'Sm'
        IF (AtomicNumber(J) == 63) AtomName(J) = 'Eu'
        IF (AtomicNumber(J) == 64) AtomName(J) = 'Gd'
        IF (AtomicNumber(J) == 65) AtomName(J) = 'Tb'
        IF (AtomicNumber(J) == 66) AtomName(J) = 'Dy'
        IF (AtomicNumber(J) == 67) AtomName(J) = 'Ho'
        IF (AtomicNumber(J) == 68) AtomName(J) = 'Er'
        IF (AtomicNumber(J) == 69) AtomName(J) = 'Tm'
        IF (AtomicNumber(J) == 70) AtomName(J) = 'Yb'
        IF (AtomicNumber(J) == 71) AtomName(J) = 'Lu'
        IF (AtomicNumber(J) == 72) AtomName(J) = 'Hf'
        IF (AtomicNumber(J) == 73) AtomName(J) = 'Ta'
        IF (AtomicNumber(J) == 74) AtomName(J) = 'W'
        IF (AtomicNumber(J) == 75) AtomName(J) = 'Re'
        IF (AtomicNumber(J) == 76) AtomName(J) = 'Os'
        IF (AtomicNumber(J) == 77) AtomName(J) = 'Ir'
        IF (AtomicNumber(J) == 78) AtomName(J) = 'Pt'
        IF (AtomicNumber(J) == 79) AtomName(J) = 'Au'
        IF (AtomicNumber(J) == 80) AtomName(J) = 'Hg'
        IF (AtomicNumber(J) == 81) AtomName(J) = 'Tl'
        IF (AtomicNumber(J) == 82) AtomName(J) = 'Pb'
        IF (AtomicNumber(J) == 83) AtomName(J) = 'Bi'
        IF (AtomicNumber(J) == 84) AtomName(J) = 'Po'
        IF (AtomicNumber(J) == 85) AtomName(J) = 'At'
        IF (AtomicNumber(J) == 86) AtomName(J) = 'Rn'
        IF (AtomicNumber(J) == 87) AtomName(J) = 'Fr'
        IF (AtomicNumber(J) == 88) AtomName(J) = 'Ra'
        IF (AtomicNumber(J) == 89) AtomName(J) = 'Ac'
        IF (AtomicNumber(J) == 90) AtomName(J) = 'Th'
        IF (AtomicNumber(J) == 91) AtomName(J) = 'Pa'
        IF (AtomicNumber(J) == 92) AtomName(J) = 'U'
        IF (AtomicNumber(J) == 93) AtomName(J) = 'Np'
        IF (AtomicNumber(J) == 94) AtomName(J) = 'Pu'
        IF (AtomicNumber(J) == 95) AtomName(J) = 'Am'
        IF (AtomicNumber(J) == 96) AtomName(J) = 'Cm'
        IF (AtomicNumber(J) == 97) AtomName(J) = 'Bk'
        IF (AtomicNumber(J) == 98) AtomName(J) = 'Cf'
        IF (AtomicNumber(J) == 99) AtomName(J) = 'Es'
        IF (AtomicNumber(J) == 100) AtomName(J) = 'Fm'
        IF (AtomicNumber(J) == 101) AtomName(J) = 'Md'
        IF (AtomicNumber(J) == 102) AtomName(J) = 'No'
        IF (AtomicNumber(J) == 103) AtomName(J) = 'Lr'
        IF (AtomicNumber(J) == 104) AtomName(J) = 'Rf'
        IF (AtomicNumber(J) == 105) AtomName(J) = 'Db'
        IF (AtomicNumber(J) == 106) AtomName(J) = 'Sg'
        IF (AtomicNumber(J) == 107) AtomName(J) = 'Bh'
        IF (AtomicNumber(J) == 108) AtomName(J) = 'Hs'
        IF (AtomicNumber(J) == 109) AtomName(J) = 'Mt'
        IF (AtomicNumber(J) == 110) AtomName(J) = 'Ds'
        IF (AtomicNumber(J) == 111) AtomName(J) = 'Rg'
        IF (AtomicNumber(J) == 112) AtomName(J) = 'Cn'
        IF (AtomicNumber(J) == 113) AtomName(J) = 'Nh'
        IF (AtomicNumber(J) == 114) AtomName(J) = 'Fl'
        IF (AtomicNumber(J) == 115) AtomName(J) = 'Mc'
        IF (AtomicNumber(J) == 116) AtomName(J) = 'Lv'
        IF (AtomicNumber(J) == 117) AtomName(J) = 'Ts'
        IF (AtomicNumber(J) == 118) AtomName(J) = 'Og'
        
        END DO
    
ELSEIF (FileType == 'Molcas') THEN
    
    ALLOCATE(AtomName(500),XYZ(500,3))
    
        AC = 1
    
        DO WHILE (.true.)
            READ(33,'(A)', END=99) Line
            Line = ADJUSTL(Line)
            IF (Line == 'Basis Set') THEN
                READ(33,*)
                DO WHILE (.true.)
                    READ(33,'(A)') Line
                    Line = ADJUSTL(Line)
                    IF (Line == 'End of Basis Set') THEN
                        EXIT
                    ELSE
                        READ(Line,*) AtomName(AC), XYZ(AC,1), XYZ(AC,2), XYZ(AC,3)
                        AC = AC + 1
                    END IF
                END DO
            END IF
        END DO
        99 CONTINUE
        NAtoms = AC - 1

END IF

!Write output file


IF (SOCounter /= 0) THEN
    OutputFile = TRIM(ADJUSTL(Cut_To_Char(InputFile,'.')))//'.xyz'
    OPEN(66, FILE = TRIM(ADJUSTL(OutputFile)), STATUS = 'UNKNOWN')
    
    
    WRITE(66,'(I0)') NAtoms
    WRITE(66,'(I0)') 
    
    DO J = 1,NAtoms
        WRITE(66,'(A5, A2, F9.6, A2, F9.6, A2, F9.6)') AtomName(J), '  ', XYZ(J,1), '  ', XYZ(J,2), '  ', XYZ(J,3)
    END DO
    
    CLOSE(66)
ELSE
    write(6,*) ' No optimised coordinates found!'
END IF
CONTAINS


FUNCTION Cut_To_Char (s,cutchar)  RESULT(outs)
    !Truncates a string at a given character
    !e.g Cut_To_Char(Hello_World,_) = Hello

CHARACTER(*)        :: s,cutchar
CHARACTER(LEN(s)+100) :: outs
INTEGER             :: i

outs = s 
DO WHILE(.TRUE.)
   i = INDEX(outs,cutchar) 
   IF (i == 0) EXIT
   outs = outs(1:i-1)
END DO
END FUNCTION Cut_To_Char

FUNCTION Cut_Past_Char (s,cutchar)  RESULT(outs)
    !Trims a string up to a given character
    !e.g Cut_Past_Char(Hello_World,_) = Hello

CHARACTER(*)        :: s,cutchar
CHARACTER(LEN(s)+100) :: outs
INTEGER             :: i

outs = s 
DO WHILE(.TRUE.)
   i = INDEX(outs,cutchar) 
   IF (i == 0) EXIT
   outs = outs(i+1:LEN(ADJUSTL(TRIM(s))))
END DO
END FUNCTION Cut_Past_Char

END PROGRAM coordinate_extractor

