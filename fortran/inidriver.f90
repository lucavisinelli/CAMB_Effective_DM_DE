    !     Code for Anisotropies in the Microwave Background
    !     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
    !     This is a sample driver routine that reads
    !     in one set of parameters and produdes the corresponding output.

    program driver
        use CAMB
        use results
    implicit none
    type (CAMBdata) :: CData
    character(len=:), allocatable :: InputFile
    integer :: i
    real(dl) :: zz(1000)

    InputFile = ''

    if (GetParamCount() /= 0)  InputFile = GetParam(1)
    if (InputFile == '') error stop 'No parameter input file'

    call CAMB_CommandLineRun(InputFile)

    end program driver

