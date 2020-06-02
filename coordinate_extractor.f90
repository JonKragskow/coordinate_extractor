PROGRAM coordinate_extractor
IMPLICIT NONE

CHARACTER(LEN = 500)              :: input_file, output_file, file_type
CHARACTER(LEN = 20)               :: result
CHARACTER(LEN = 100), ALLOCATABLE :: atom_symbols(:)
INTEGER                           :: num_atoms
LOGICAL                           :: geom_found
REAL(KIND = 8), ALLOCATABLE       :: XYZ(:,:)


!Read user input
CALL READ_USER_INPUT(input_file, file_type, output_file)

!Read input file
IF (file_type == 'Gaussianout' .OR. file_type == 'Gaussianin' ) THEN
    CALL READ_GAUSSIAN(input_file, file_type, num_atoms, XYZ, atom_symbols, result, geom_found)
ELSE IF (file_type == 'Molcas') THEN
    CALL READ_MOLCAS(input_file, file_type, num_atoms, XYZ, atom_symbols, result, geom_found)
END IF

!Write xyz file
IF (geom_found .EQV. .TRUE.) THEN
    CALL WRITE_XYZ(num_atoms, XYZ, atom_symbols, output_file, result)
END IF


CONTAINS

SUBROUTINE READ_USER_INPUT(input_file, file_type, output_file)
    IMPLICIT NONE
    CHARACTER(LEN = *), INTENT(OUT) :: input_file, file_type, output_file
    CHARACTER(LEN = 10)             :: file_extension
    LOGICAL                         :: file_exists

    CALL GET_COMMAND_ARGUMENT(1,input_file)
    
    IF (input_file == '-h' .OR. input_file == '') THEN
        WRITE(6,'(A)') 'coordinate_extractor <file>'
        WRITE(6,'(A)')
        WRITE(6,'(A)') 'Extracts coordinates from a Gaussian log or molcas .out/.in file'
        WRITE(6,'(A)') 'file             : CHARACTER          output file of calculation                    e.g  dy.log'
        STOP
    END IF

    WRITE(file_extension,'(A)') ADJUSTL(TRIM(LOWERCASE(Cut_Past_Char(input_file,'.'))))

    IF (file_extension(1:2) == 'lo') THEN
        file_type = 'Gaussianout'
        WRITE(output_file,'(A, A)') ADJUSTL(TRIM(Cut_To_Char(input_file,'.'))), '_opt.xyz'
    ELSE IF (file_extension(1:2) == 'co') THEN
        file_type = 'Gaussianin'
        WRITE(output_file,'(A, A)') ADJUSTL(TRIM(Cut_To_Char(input_file,'.'))), '.xyz'
    ELSE IF (file_extension(1:2) == 'ou' .OR. file_extension(1:2) == 'in') THEN
        file_type = 'Molcas'
        WRITE(output_file,'(A, A)') ADJUSTL(TRIM(Cut_To_Char(input_file,'.'))), '.xyz'
    ELSE
        WRITE(6,*) 'Unknown input file type'
        STOP
    END IF

    !Check file exists
    INQUIRE(FILE=input_file, EXIST=file_exists)
    IF (file_exists .EQV. .FALSE.) THEN
        WRITE(6,'(3A)') 'Specified file', TRIM(ADJUSTL(input_file)),' does not exist'
        STOP
    END IF

END SUBROUTINE READ_USER_INPUT

SUBROUTINE READ_GAUSSIAN(input_file, file_type, num_atoms, XYZ, atom_symbols, result, geom_found)
    IMPLICIT NONE
    CHARACTER(LEN = *), INTENT(IN)                 :: input_file, file_type
    CHARACTER(LEN = 500)                           :: line, CDummy
    CHARACTER(LEN = 20), INTENT(OUT)               :: result
    CHARACTER(LEN = 100), ALLOCATABLE, INTENT(OUT) :: atom_symbols(:)
    INTEGER                                        :: reason, count_SO, doubledash, Idummy, J
    INTEGER, INTENT(OUT)                           :: num_atoms
    LOGICAL                                        :: allocated, zmat
    LOGICAL, INTENT(OUT)                           :: geom_found
    INTEGER, ALLOCATABLE                           :: atomic_numbers(:)
    REAL(KIND = 8), ALLOCATABLE, INTENT(OUT)       :: XYZ(:,:)

    count_SO    = 0
    geom_found  = .FALSE.
    allocated = .FALSE.
    doubledash  = 0
    zmat = .FALSE.

    OPEN(33, FILE = TRIM(ADJUSTL(input_file)), STATUS = 'OLD')
    
    IF (file_type == 'Gaussianout') THEN

        DO
            !Read in a line from the file
            READ(33,'(A)',IOSTAT=reason) line
            line = ADJUSTL(line)
        
            !Check iostat output for file
            IF (Reason > 0)  THEN
                WRITE(6,'(A)') 'UNKNOWN ERROR'
                WRITE(6,'(A)') 'ABORTING.....'
            ELSE IF (Reason < 0) THEN
                EXIT
            ELSE

                !Get number of atoms and allocate arrays accordingly
                IF (line(1:7) == 'NAtoms=') THEN
                    READ(line,*) CDummy, num_atoms
                    IF (allocated .EQV. .FALSE.) THEN
                        ALLOCATE(atomic_numbers(num_atoms), XYZ(num_atoms,3), atom_symbols(num_atoms))
                        atomic_numbers = 0
                        XYZ = 0.0_8
                        allocated = .TRUE.
                    END IF
                END IF
    
                !Read in coordinates and labels
                IF (allocated .EQV. .TRUE.) THEN
        
                    !Read in coordinates from checkpoint geometry
                    IF (line(1:35) == 'Structure from the checkpoint file:') THEN
                        DO WHILE (doubledash /= 2)
                            READ(33,*) line
                            IF (line(1:5) == '-----') doubledash = doubledash + 1
                        END DO
                        READ(33,*)
                        READ(33,*) line
                        IF (line(1:8) == 'Z-Matrix') THEN
                            zmat = .TRUE.
                        END IF
                        
                        IF (zmat .EQV. .FALSE.) THEN
                            DO J = 1, num_atoms
                                READ(33,*) atom_symbols(J), IDummy, XYZ(J,1), XYZ(J,2), XYZ(J,3)
                            END DO
                            geom_found = .TRUE.
                        END IF
                    END IF
                    
                    !Read in coordinates from standard orientation geometry
                    IF (line(1:21) == 'Standard orientation:') THEN
                            READ(33,*)
                            READ(33,*)
                            READ(33,*)
                            READ(33,*)
                            DO J = 1,num_atoms
                                READ(33,*) Idummy, atomic_numbers(J), Idummy, XYZ(J,1), XYZ(J,2), XYZ(J,3)
                            END DO
                        !Keep track of number of standard orientation sections
                        count_SO = count_SO + 1
                        geom_found = .TRUE.
                    END IF
                END IF
            END IF
        END DO 
            
        !Close input file
        CLOSE(33)
            
        !Replace Atomic Numbers with Letters in new array
        DO J = 1, num_atoms
        
            IF (atomic_numbers(J) == 1) atom_symbols(J) = 'H'
            IF (atomic_numbers(J) == 2) atom_symbols(J) = 'He'
            IF (atomic_numbers(J) == 3) atom_symbols(J) = 'Li'
            IF (atomic_numbers(J) == 4) atom_symbols(J) = 'Be'
            IF (atomic_numbers(J) == 5) atom_symbols(J) = 'B'
            IF (atomic_numbers(J) == 6) atom_symbols(J) = 'C'
            IF (atomic_numbers(J) == 7) atom_symbols(J) = 'N'
            IF (atomic_numbers(J) == 8) atom_symbols(J) = 'O'
            IF (atomic_numbers(J) == 9) atom_symbols(J) = 'F'
            IF (atomic_numbers(J) == 10) atom_symbols(J) = 'Ne'
            IF (atomic_numbers(J) == 11) atom_symbols(J) = 'Na'
            IF (atomic_numbers(J) == 12) atom_symbols(J) = 'Mg'
            IF (atomic_numbers(J) == 13) atom_symbols(J) = 'Al'
            IF (atomic_numbers(J) == 14) atom_symbols(J) = 'Si'
            IF (atomic_numbers(J) == 15) atom_symbols(J) = 'P'
            IF (atomic_numbers(J) == 16) atom_symbols(J) = 'S'
            IF (atomic_numbers(J) == 17) atom_symbols(J) = 'Cl'
            IF (atomic_numbers(J) == 18) atom_symbols(J) = 'Ar'
            IF (atomic_numbers(J) == 19) atom_symbols(J) = 'K'
            IF (atomic_numbers(J) == 20) atom_symbols(J) = 'Ca'
            IF (atomic_numbers(J) == 21) atom_symbols(J) = 'Sc'
            IF (atomic_numbers(J) == 22) atom_symbols(J) = 'Ti'
            IF (atomic_numbers(J) == 23) atom_symbols(J) = 'V'
            IF (atomic_numbers(J) == 24) atom_symbols(J) = 'Cr'
            IF (atomic_numbers(J) == 25) atom_symbols(J) = 'Mn'
            IF (atomic_numbers(J) == 26) atom_symbols(J) = 'Fe'
            IF (atomic_numbers(J) == 27) atom_symbols(J) = 'Co'
            IF (atomic_numbers(J) == 28) atom_symbols(J) = 'Ni'
            IF (atomic_numbers(J) == 29) atom_symbols(J) = 'Cu'
            IF (atomic_numbers(J) == 30) atom_symbols(J) = 'Zn'
            IF (atomic_numbers(J) == 31) atom_symbols(J) = 'Ga'
            IF (atomic_numbers(J) == 32) atom_symbols(J) = 'Ge'
            IF (atomic_numbers(J) == 33) atom_symbols(J) = 'As'
            IF (atomic_numbers(J) == 34) atom_symbols(J) = 'Se'
            IF (atomic_numbers(J) == 35) atom_symbols(J) = 'Br'
            IF (atomic_numbers(J) == 36) atom_symbols(J) = 'Kr'
            IF (atomic_numbers(J) == 37) atom_symbols(J) = 'Rb'
            IF (atomic_numbers(J) == 38) atom_symbols(J) = 'Sr'
            IF (atomic_numbers(J) == 39) atom_symbols(J) = 'Y'
            IF (atomic_numbers(J) == 40) atom_symbols(J) = 'Zr'
            IF (atomic_numbers(J) == 41) atom_symbols(J) = 'Nb'
            IF (atomic_numbers(J) == 42) atom_symbols(J) = 'Mo'
            IF (atomic_numbers(J) == 43) atom_symbols(J) = 'Tc'
            IF (atomic_numbers(J) == 44) atom_symbols(J) = 'Ru'
            IF (atomic_numbers(J) == 45) atom_symbols(J) = 'Rh'
            IF (atomic_numbers(J) == 46) atom_symbols(J) = 'Pd'
            IF (atomic_numbers(J) == 47) atom_symbols(J) = 'Ag'
            IF (atomic_numbers(J) == 48) atom_symbols(J) = 'Cd'
            IF (atomic_numbers(J) == 49) atom_symbols(J) = 'In'
            IF (atomic_numbers(J) == 50) atom_symbols(J) = 'Sn'
            IF (atomic_numbers(J) == 51) atom_symbols(J) = 'Sb'
            IF (atomic_numbers(J) == 52) atom_symbols(J) = 'Te'
            IF (atomic_numbers(J) == 53) atom_symbols(J) = 'I'
            IF (atomic_numbers(J) == 54) atom_symbols(J) = 'Xe'
            IF (atomic_numbers(J) == 55) atom_symbols(J) = 'Cs'
            IF (atomic_numbers(J) == 56) atom_symbols(J) = 'Ba'
            IF (atomic_numbers(J) == 57) atom_symbols(J) = 'La'
            IF (atomic_numbers(J) == 58) atom_symbols(J) = 'Ce'
            IF (atomic_numbers(J) == 59) atom_symbols(J) = 'Pr'
            IF (atomic_numbers(J) == 60) atom_symbols(J) = 'Nd'
            IF (atomic_numbers(J) == 61) atom_symbols(J) = 'Pm'
            IF (atomic_numbers(J) == 62) atom_symbols(J) = 'Sm'
            IF (atomic_numbers(J) == 63) atom_symbols(J) = 'Eu'
            IF (atomic_numbers(J) == 64) atom_symbols(J) = 'Gd'
            IF (atomic_numbers(J) == 65) atom_symbols(J) = 'Tb'
            IF (atomic_numbers(J) == 66) atom_symbols(J) = 'Dy'
            IF (atomic_numbers(J) == 67) atom_symbols(J) = 'Ho'
            IF (atomic_numbers(J) == 68) atom_symbols(J) = 'Er'
            IF (atomic_numbers(J) == 69) atom_symbols(J) = 'Tm'
            IF (atomic_numbers(J) == 70) atom_symbols(J) = 'Yb'
            IF (atomic_numbers(J) == 71) atom_symbols(J) = 'Lu'
            IF (atomic_numbers(J) == 72) atom_symbols(J) = 'Hf'
            IF (atomic_numbers(J) == 73) atom_symbols(J) = 'Ta'
            IF (atomic_numbers(J) == 74) atom_symbols(J) = 'W'
            IF (atomic_numbers(J) == 75) atom_symbols(J) = 'Re'
            IF (atomic_numbers(J) == 76) atom_symbols(J) = 'Os'
            IF (atomic_numbers(J) == 77) atom_symbols(J) = 'Ir'
            IF (atomic_numbers(J) == 78) atom_symbols(J) = 'Pt'
            IF (atomic_numbers(J) == 79) atom_symbols(J) = 'Au'
            IF (atomic_numbers(J) == 80) atom_symbols(J) = 'Hg'
            IF (atomic_numbers(J) == 81) atom_symbols(J) = 'Tl'
            IF (atomic_numbers(J) == 82) atom_symbols(J) = 'Pb'
            IF (atomic_numbers(J) == 83) atom_symbols(J) = 'Bi'
            IF (atomic_numbers(J) == 84) atom_symbols(J) = 'Po'
            IF (atomic_numbers(J) == 85) atom_symbols(J) = 'At'
            IF (atomic_numbers(J) == 86) atom_symbols(J) = 'Rn'
            IF (atomic_numbers(J) == 87) atom_symbols(J) = 'Fr'
            IF (atomic_numbers(J) == 88) atom_symbols(J) = 'Ra'
            IF (atomic_numbers(J) == 89) atom_symbols(J) = 'Ac'
            IF (atomic_numbers(J) == 90) atom_symbols(J) = 'Th'
            IF (atomic_numbers(J) == 91) atom_symbols(J) = 'Pa'
            IF (atomic_numbers(J) == 92) atom_symbols(J) = 'U'
            IF (atomic_numbers(J) == 93) atom_symbols(J) = 'Np'
            IF (atomic_numbers(J) == 94) atom_symbols(J) = 'Pu'
            IF (atomic_numbers(J) == 95) atom_symbols(J) = 'Am'
            IF (atomic_numbers(J) == 96) atom_symbols(J) = 'Cm'
            IF (atomic_numbers(J) == 97) atom_symbols(J) = 'Bk'
            IF (atomic_numbers(J) == 98) atom_symbols(J) = 'Cf'
            IF (atomic_numbers(J) == 99) atom_symbols(J) = 'Es'
            IF (atomic_numbers(J) == 100) atom_symbols(J) = 'Fm'
            IF (atomic_numbers(J) == 101) atom_symbols(J) = 'Md'
            IF (atomic_numbers(J) == 102) atom_symbols(J) = 'No'
            IF (atomic_numbers(J) == 103) atom_symbols(J) = 'Lr'
            IF (atomic_numbers(J) == 104) atom_symbols(J) = 'Rf'
            IF (atomic_numbers(J) == 105) atom_symbols(J) = 'Db'
            IF (atomic_numbers(J) == 106) atom_symbols(J) = 'Sg'
            IF (atomic_numbers(J) == 107) atom_symbols(J) = 'Bh'
            IF (atomic_numbers(J) == 108) atom_symbols(J) = 'Hs'
            IF (atomic_numbers(J) == 109) atom_symbols(J) = 'Mt'
            IF (atomic_numbers(J) == 110) atom_symbols(J) = 'Ds'
            IF (atomic_numbers(J) == 111) atom_symbols(J) = 'Rg'
            IF (atomic_numbers(J) == 112) atom_symbols(J) = 'Cn'
            IF (atomic_numbers(J) == 113) atom_symbols(J) = 'Nh'
            IF (atomic_numbers(J) == 114) atom_symbols(J) = 'Fl'
            IF (atomic_numbers(J) == 115) atom_symbols(J) = 'Mc'
            IF (atomic_numbers(J) == 116) atom_symbols(J) = 'Lv'
            IF (atomic_numbers(J) == 117) atom_symbols(J) = 'Ts'
            IF (atomic_numbers(J) == 118) atom_symbols(J) = 'Og'
            
        END DO
        IF (geom_found .EQV. .TRUE. .AND. count_SO > 0 ) WRITE(result,'(A)') 'optimisation'
        IF (geom_found .EQV. .TRUE. .AND. count_SO == 0) WRITE(result,'(A)') 'checkpoint'

    ELSE IF (file_type == 'Gaussianin') THEN
    
        WRITE(6,'(A)') '----  Currently unsupported   ----'
    
        STOP

    END IF
    
END SUBROUTINE READ_GAUSSIAN

SUBROUTINE READ_MOLCAS(input_file, file_type, num_atoms, XYZ, atom_symbols, result, geom_found)
    IMPLICIT NONE
    CHARACTER(LEN = *), INTENT(IN)                 :: input_file, file_type
    CHARACTER(LEN = 500)                           :: line
    CHARACTER(LEN = 20), INTENT(OUT)               :: result
    CHARACTER(LEN = 100), ALLOCATABLE, INTENT(OUT) :: atom_symbols(:)
    INTEGER                                        :: reason, num_atoms, count
    LOGICAL, INTENT(OUT)                           :: geom_found
    REAL(KIND = 8), ALLOCATABLE, INTENT(OUT)       :: XYZ(:,:)

    geom_found  = .FALSE.
    count  = 0


    OPEN(33, FILE = TRIM(ADJUSTL(input_file)), STATUS = 'OLD')
    
        ALLOCATE(atom_symbols(500),XYZ(500,3))
        
        count = 1
    
        DO
            !Read in a line from the file
            READ(33,'(A)',IOSTAT=reason) line
            line = ADJUSTL(line)
        
            !Check iostat output for file
            IF (Reason > 0)  THEN
                WRITE(6,'(A)') 'UNKNOWN ERROR'
                WRITE(6,'(A)') 'ABORTING.....'
            ELSE IF (Reason < 0) THEN
                EXIT
            ELSE
                IF (line == 'Basis Set') THEN
                    READ(33,*)
                    DO WHILE (.true.)
                        READ(33,'(A)') line
                        line = ADJUSTL(line)
                        IF (line == 'End of Basis Set') THEN
                            EXIT
                        ELSE
                            READ(line,*) atom_symbols(count), XYZ(count,1), XYZ(count,2), XYZ(count,3)
                            count = count + 1
                        END IF
                    END DO
                END IF
            END IF
        END DO
    CLOSE(33)
    
    num_atoms = count - 1
    IF (num_atoms > 0) geom_found = .TRUE.

    WRITE(result,'(A)') 'input'

END SUBROUTINE READ_MOLCAS

SUBROUTINE WRITE_XYZ(num_atoms, XYZ, atom_symbols, output_file, result)
    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: num_atoms
    INTEGER                        :: J, num_stars
    REAL(KIND=8), INTENT(IN)       :: XYZ(:,:)
    CHARACTER(LEN = *), INTENT(IN) :: atom_symbols(:), output_file, result
    CHARACTER(LEN = 1000)          :: output_char, stars, FMT

    OPEN(66, FILE = TRIM(ADJUSTL(output_file)), STATUS = 'UNKNOWN')
        
        WRITE(66,'(I0)') num_atoms
        WRITE(66,*) 
        
        DO J = 1,num_atoms
            WRITE(66,'(A, A, F9.6, A, F9.6, A, F9.6)') ADJUSTL(TRIM(atom_symbols(J))), '  ', XYZ(J,1), '  ', XYZ(J,2), '  ', XYZ(J,3)
        END DO
    
    CLOSE(66)
    

    !Output info on xyz file to terminal

    WRITE(output_char,'(A, A, A, A, A)') '---- Coordinates from ',TRIM(ADJUSTL(result)),' written to ', TRIM(ADJUSTL(output_file)),' ----'
    num_stars = LEN(TRIM(ADJUSTL(output_char)))
    WRITE(FMT,'(A,I0, A)') '(',num_stars,'A)'
    WRITE(6,ADJUSTL(TRIM(FMT))) ('*',J = 1, num_stars)
    WRITE(6,'(A)') TRIM(ADJUSTL(output_char))
    WRITE(6,ADJUSTL(TRIM(FMT))) ('*',J = 1, num_stars)

END SUBROUTINE WRITE_XYZ

FUNCTION LOWERCASE(in_str) RESULT(out_str)
    !Converts str to lowercase characters
    IMPLICIT NONE
    INTEGER                           :: i,gap
    CHARACTER(LEN = *), INTENT(IN) :: in_str
    CHARACTER(LEN = LEN(in_str))      :: out_str

    out_str = in_str
    
    gap = ICHAR('a')-ICHAR('A')
    IF(len(out_str) > 0) THEN
        DO i=1,len(out_str)
            IF(out_str(i:i) <= 'Z') THEN
                IF(out_str(i:i) >= 'A') out_str(i:i) = CHAR(ICHAR(out_str(i:i)) + gap)
            END IF
        END DO
END IF
END FUNCTION LOWERCASE

FUNCTION Cut_To_Char (s,cutchar)  RESULT(outs)
    !Truncates a string at a given character
    !e.g Cut_To_Char('Hello_World','_') = Hello

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
    !Trims a string past to a given character
    !e.g Cut_Past_Char('Hello_World','_') = World

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

