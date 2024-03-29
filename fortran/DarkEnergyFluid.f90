
    module DarkEnergyFluid
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none

    type, extends(TDarkEnergyEqnOfState) :: TDarkEnergyFluid
        !comoving sound speed is always exactly 1 for quintessence
        !(otherwise assumed constant, though this is almost certainly unrealistic)
    contains
    procedure :: ReadParams => TDarkEnergyFluid_ReadParams
    procedure, nopass :: PythonClass => TDarkEnergyFluid_PythonClass
    procedure, nopass :: SelfPointer => TDarkEnergyFluid_SelfPointer
    procedure :: Init =>TDarkEnergyFluid_Init
    procedure :: PerturbedStressEnergy => TDarkEnergyFluid_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TDarkEnergyFluid_PerturbationEvolve
    end type TDarkEnergyFluid

    !Example implementation of fluid model using specific analytic form
    !(approximate effective axion fluid model from arXiv:1806.10608, with c_s^2=1 if n=infinity (w_n=1))
    !This is an example, it's not supposed to be a rigorous model!  (not very well tested)
    type, extends(TDarkEnergyModel) :: TAxionEffectiveFluid
        real(dl) :: w_n = 1._dl !Effective equation of state when oscillating
        real(dl) :: Om = 0._dl !Omega of the early DE component today (assumed to be negligible compared to omega_lambda)
        real(dl) :: a_c  !transition scale factor
        real(dl) :: theta_i = const_pi/2 !Initial value
        real(dl), private :: pow, omL, acpow, freq, n !cached internally
    contains
    procedure :: ReadParams =>  TAxionEffectiveFluid_ReadParams
    procedure, nopass :: PythonClass => TAxionEffectiveFluid_PythonClass
    procedure, nopass :: SelfPointer => TAxionEffectiveFluid_SelfPointer
    procedure :: Init => TAxionEffectiveFluid_Init
    procedure :: w_de => TAxionEffectiveFluid_w_de
    procedure :: grho_de => TAxionEffectiveFluid_grho_de
    procedure :: PerturbedStressEnergy => TAxionEffectiveFluid_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TAxionEffectiveFluid_PerturbationEvolve
    end type TAxionEffectiveFluid

    type, extends(TDarkEnergyFluid) :: TDMDEInteraction
    !    real(dl) :: w_width = 0.1_dl !may not be used
    !    real(dl) :: a_dec   = 0.1_dl !may not be used
    contains
    procedure :: ReadParams =>  TDMDEInteraction_ReadParams
    procedure, nopass :: PythonClass => TDMDEInteraction_PythonClass
    procedure, nopass :: SelfPointer => TDMDEInteraction_SelfPointer
    procedure :: Init => TDMDEInteraction_Init
    procedure :: w_de => TDMDEInteraction_w_de
    procedure :: grho_de => TDMDEInteraction_grho_de
    procedure :: PerturbationEvolve => TDMDEInteraction_PerturbationEvolve
    procedure :: PerturbedStressEnergy => TDMDEInteraction_PerturbedStressEnergy
    end type TDMDEInteraction

    contains


    subroutine TDarkEnergyFluid_ReadParams(this, Ini)
    use IniObjects
    class(TDarkEnergyFluid) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyEqnOfState%ReadParams(Ini)
    this%cs2_lam = Ini%Read_Double('cs2_lam', 1.d0)
    end subroutine TDarkEnergyFluid_ReadParams


    function TDarkEnergyFluid_PythonClass()
    character(LEN=:), allocatable :: TDarkEnergyFluid_PythonClass

    TDarkEnergyFluid_PythonClass = 'DarkEnergyFluid'

    end function TDarkEnergyFluid_PythonClass

    subroutine TDarkEnergyFluid_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TDarkEnergyFluid), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TDarkEnergyFluid_SelfPointer

    subroutine TDarkEnergyFluid_Init(this, State)
    use classes
    class(TDarkEnergyFluid), intent(inout) :: this
    class(TCAMBdata), intent(in) :: State

    call this%TDarkEnergyEqnOfState%Init(State)

    if (this%is_cosmological_constant) then
        this%num_perturb_equations = 0
    else
        if (this%use_tabulated_w) then
            if (any(this%equation_of_state%F<-1)) &
                error stop 'Fluid dark energy model does not allow w crossing -1'
        elseif (this%wa/=0 .and. &
            ((1+this%w_lam < -1.e-6_dl) .or. 1+this%w_lam + this%wa < -1.e-6_dl)) then
            error stop 'Fluid dark energy model does not allow w crossing -1'
        end if
        this%num_perturb_equations = 2
    end if

    end subroutine TDarkEnergyFluid_Init


    subroutine TDarkEnergyFluid_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDarkEnergyFluid), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    if (this%no_perturbations) then
        dgrhoe=0
        dgqe=0
    else
        dgrhoe = ay(w_ix) * grhov_t
        dgqe = ay(w_ix + 1) * grhov_t * (1 + w)
    end if
    end subroutine TDarkEnergyFluid_PerturbedStressEnergy


    subroutine TDarkEnergyFluid_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    class(TDarkEnergyFluid), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) Hv3_over_k, loga
    real(dl) :: weff, adecay, width

    Hv3_over_k =  3*adotoa* y(w_ix + 1) / k
    !density perturbation
    ayprime(w_ix) = -3 * adotoa * (this%cs2_lam - w) *  (y(w_ix) + (1 + w) * Hv3_over_k) &
        -  (1 + w) * k * y(w_ix + 1) - (1 + w) * k * z
    if (this%use_tabulated_w) then
        !account for derivatives of w
        loga = log(a)
        if (loga > this%equation_of_state%Xmin_interp .and. loga < this%equation_of_state%Xmax_interp) then
            ayprime(w_ix) = ayprime(w_ix) - adotoa*this%equation_of_state%Derivative(loga)* Hv3_over_k
        end if
    elseif (this%wa/=0) then
        ayprime(w_ix) = ayprime(w_ix) + Hv3_over_k*this%wa*adotoa*a
    end if
    !velocity
    if (abs(w+1) > 1e-6) then
        ayprime(w_ix + 1) = -adotoa * (1 - 3 * this%cs2_lam) * y(w_ix + 1) + &
            k * this%cs2_lam * y(w_ix) / (1 + w)
    else
        ayprime(w_ix + 1) = 0
    end if

    end subroutine TDarkEnergyFluid_PerturbationEvolve



    subroutine TAxionEffectiveFluid_ReadParams(this, Ini)
    use IniObjects
    class(TAxionEffectiveFluid) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyModel%ReadParams(Ini)
    this%w_n = Ini%Read_Double('AxionEffectiveFluid_w_n')
    this%om  = Ini%Read_Double('AxionEffectiveFluid_om')
    this%a_c = Ini%Read_Double('AxionEffectiveFluid_a_c')
    call Ini%Read('AxionEffectiveFluid_theta_i', this%theta_i)

    end subroutine TAxionEffectiveFluid_ReadParams


    function TAxionEffectiveFluid_PythonClass()
    character(LEN=:), allocatable :: TAxionEffectiveFluid_PythonClass

    TAxionEffectiveFluid_PythonClass = 'AxionEffectiveFluid'
    end function TAxionEffectiveFluid_PythonClass

    subroutine TAxionEffectiveFluid_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TAxionEffectiveFluid), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TAxionEffectiveFluid_SelfPointer

    subroutine TAxionEffectiveFluid_Init(this, State)
    use classes
    class(TAxionEffectiveFluid), intent(inout) :: this
    class(TCAMBdata), intent(in) :: State
    real(dl) :: grho_rad, F, p, mu, xc, n

    select type(State)
    class is (CAMBdata)
        this%is_cosmological_constant = this%om==0
        this%pow = 3*(1+this%w_n)
        this%omL = State%Omega_de - this%om !Omega_de is total dark energy density today
        this%acpow = this%a_c**this%pow
        this%num_perturb_equations = 2
        if (this%w_n < 0.9999) then
            ! n <> infinity
            !get (very) approximate result for sound speed parameter; arXiv:1806.10608  Eq 30 (but mu may not exactly agree with what they used)
            n = nint((1+this%w_n)/(1-this%w_n))
            !Assume radiation domination, standard neutrino model; H0 factors cancel
            grho_rad = (kappa/c**2*4*sigma_boltz/c**3*State%CP%tcmb**4*Mpc**2*(1+3.046*7._dl/8*(4._dl/11)**(4._dl/3)))
            xc = this%a_c**2/2/sqrt(grho_rad/3)
            F=7./8
            p=1./2
            mu = 1/xc*(1-cos(this%theta_i))**((1-n)/2.)*sqrt((1-F)*(6*p+2)*this%theta_i/n/sin(this%theta_i))
            this%freq =  mu*(1-cos(this%theta_i))**((n-1)/2.)* &
                sqrt(const_pi)*Gamma((n+1)/(2.*n))/Gamma(1+0.5/n)*2.**(-(n**2+1)/(2.*n))*3.**((1./n-1)/2)*this%a_c**(-6./(n+1)+3) &
                *( this%a_c**(6*n/(n+1.))+1)**(0.5*(1./n-1))
            this%n = n
        end if
    end select

    end subroutine TAxionEffectiveFluid_Init


    function TAxionEffectiveFluid_w_de(this, a)
    class(TAxionEffectiveFluid) :: this
    real(dl) :: TAxionEffectiveFluid_w_de
    real(dl), intent(IN) :: a
    real(dl) :: rho, apow, acpow

    apow = a**this%pow
    acpow = this%acpow
    rho = this%omL+ this%om*(1+acpow)/(apow+acpow)
    TAxionEffectiveFluid_w_de = this%om*(1+acpow)/(apow+acpow)**2*(1+this%w_n)*apow/rho - 1

    end function TAxionEffectiveFluid_w_de

    function TAxionEffectiveFluid_grho_de(this, a) result(grho_de)  !relative density (8 pi G a^4 rho_de /grhov)
      class(TAxionEffectiveFluid) :: this
      real(dl) :: grho_de, apow
      real(dl), intent(IN) :: a

      if(a == 0.d0)then
        grho_de = 0.d0
      else
        apow = a**this%pow
        grho_de = (this%omL*(apow+this%acpow)+this%om*(1+this%acpow))*a**4 &
            /((apow+this%acpow)*(this%omL+this%om))
      endif
    end function TAxionEffectiveFluid_grho_de

    subroutine TAxionEffectiveFluid_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    class(TAxionEffectiveFluid), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) Hv3_over_k, deriv, apow, acpow, cs2, fac

    if (this%w_n < 0.9999) then
        fac = 2*a**(2-6*this%w_n)*this%freq**2
        cs2 = (fac*(this%n-1) + k**2)/(fac*(this%n+1) + k**2)
    else
        cs2 = 1
    end if
    apow = a**this%pow
    acpow = this%acpow
    Hv3_over_k =  3*adotoa* y(w_ix + 1) / k
    ! dw/dlog a/(1+w)
    deriv  = (acpow**2*(this%om+this%omL)+this%om*acpow-apow**2*this%omL)*this%pow &
        /((apow+acpow)*(this%omL*(apow+acpow)+this%om*(1+acpow)))
    !density perturbation
    ayprime(w_ix) = -3 * adotoa * (cs2 - w) *  (y(w_ix) + Hv3_over_k) &
        -   k * y(w_ix + 1) - (1 + w) * k * z  - adotoa*deriv* Hv3_over_k
    !(1+w)v
    ayprime(w_ix + 1) = -adotoa * (1 - 3 * cs2 - deriv) * y(w_ix + 1) + &
        k * cs2 * y(w_ix)

    end subroutine TAxionEffectiveFluid_PerturbationEvolve


    subroutine TAxionEffectiveFluid_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TAxionEffectiveFluid), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe = ay(w_ix) * grhov_t
    dgqe = ay(w_ix + 1) * grhov_t

    end subroutine TAxionEffectiveFluid_PerturbedStressEnergy

    !!!!!!!!!!
    !! TDMDEInteraction
    !!

    subroutine TDMDEInteraction_ReadParams(this, Ini)
    use IniObjects
    class(TDMDEInteraction)  :: this
    class(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyEqnOfState%ReadParams(Ini)
      this%w_width = Ini%Read_Double('w_width', 0.5d0)
      this%a_dec   = Ini%Read_Double('a_dec',   0.1d0)
      this%m_phi   = Ini%Read_Double('m_phi',   1.d-27)
    end subroutine TDMDEInteraction_ReadParams

    function TDMDEInteraction_PythonClass()
    character(LEN=:), allocatable :: TDMDEInteraction_PythonClass
    TDMDEInteraction_PythonClass = 'DMDEInteraction'
    end function TDMDEInteraction_PythonClass

    subroutine TDMDEInteraction_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TDMDEInteraction), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P
    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine TDMDEInteraction_SelfPointer

    subroutine TDMDEInteraction_Init(this, State)
    use classes
    class(TDMDEInteraction), intent(inout) :: this
    class(TCAMBdata), intent(in) :: State

        call this%TDarkEnergyEqnOfState%Init(State)

        if (this%is_cosmological_constant) then
          this%num_perturb_equations = 0
        else
          this%num_perturb_equations = 2
        end if

    end subroutine TDMDEInteraction_Init

    function TDMDEInteraction_w_de(this, a)
      class(TDMDEInteraction) :: this
      real(dl) :: TDMDEInteraction_w_de
      real(dl), intent(IN) :: a
      real(dl) :: weff, adecay, width

      weff   = this%w_lam
      adecay = this%a_dec
      width  = this%w_width
      TDMDEInteraction_w_de = 0.5*weff*(1 + TANH(LOG(a/adecay)/width))

    end function TDMDEInteraction_w_de

    function TDMDEInteraction_grho_de(this, a) result(grho_de)
    class(TDMDEInteraction) :: this
    real(dl) :: grho_de
    real(dl), intent(IN) :: a
    real(dl) :: Log0, Log1
    real(dl) :: weff, adecay, width

        weff   = this%w_lam
        adecay = this%a_dec
        width  = this%w_width

       !! We use w(a) = weff/2 (1 + Tanh[Log[a/adecay]/width])
       !! Here we compute
       !!
       !!  exp( Integrate( -3(1+w) da/a) )* a**4

       Log0 = LOG(1/adecay)/width
       Log1 = LOG(a/adecay)/width
       grho_de = a**(1 - 1.5_dl*weff)*EXP( 1.5_dl*weff*width*LOG( COSH(Log0) / COSH(Log1) ) )
    end function TDMDEInteraction_grho_de

    subroutine TDMDEInteraction_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
      class(TDMDEInteraction), intent(in) :: this
      real(dl), intent(inout) :: ayprime(:)
      real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
      integer, intent(in) :: w_ix
      real(dl), parameter :: Mpcm1 = 6.4e-30_dl 
      real(dl) Hv3_over_k, cs2, fac, wp1, wprime
      real(dl) :: weff, adecay, width, mphi

      weff   = this%w_lam
      adecay = this%a_dec
      width  = this%w_width
      mphi   = this%m_phi

      !!
      !! cs2 = k^2 / (k^2 + 4m^2 a^2)
      !! The axion mass is in units of hbar c / Mpc = 6.4*10^-30eV
      !! In the following, we assume a mass 10^-27eV
      !! so there is an extra factor 300

      if (a < adecay) then
        fac = 2*a*mphi/Mpcm1
        cs2 = k**2/(k**2 + fac**2 )
      else
        cs2 = 1
      end if
      !write(*,*) cs2
      !cs2=abs(w)

      !!
      !! We set u = (1+w)v and we solve Eqs.35-36 in the CAMB notebook
      !! We account for the derivative of w as wprime=dw/dlog a/(1+w)
      !!
      wprime = weff/(width* COSH( LOG(a/adecay)/width)**2) / (2*(1+w))

      !! Eq.35 in the CAMB notebook for a varying equation of state
      Hv3_over_k =  3*adotoa* y(w_ix + 1) / k
      ayprime(w_ix) = -3 * adotoa * (cs2 - w) * (y(w_ix) + Hv3_over_k) &
            -  k * y(w_ix + 1) - (1 + w) * k * z - adotoa*wprime*Hv3_over_k

      !! Eq.36 for u=(1+w)v
      ayprime(w_ix + 1) = -adotoa * (1 - 3 * cs2 - wprime) * y(w_ix + 1) + k * cs2 * y(w_ix)

    end subroutine TDMDEInteraction_PerturbationEvolve

    subroutine TDMDEInteraction_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TDMDEInteraction), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    dgrhoe = ay(w_ix) * grhov_t
    dgqe = ay(w_ix + 1) * grhov_t

    end subroutine TDMDEInteraction_PerturbedStressEnergy


    end module DarkEnergyFluid
