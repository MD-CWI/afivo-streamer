!> Module to set the type of model
module m_model
  use m_af_all
  use m_types

  implicit none
  private

  !> Fluid model with local field approximation
  integer, parameter :: model_lfa = 1

  !> Fluid model with local energy approximation and energy fluxes that are
  !> "5/3" times the electron flux
  integer, parameter :: model_ee53 = 2

  !> Which type of model is used
  integer, public, protected :: model_type = model_lfa

  !> Whether the model has an energy equation
  logical, public, protected :: model_has_energy_equation = .false.

  public :: model_initialize

contains

  !> Initialize the module
  subroutine model_initialize(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg
    character(len=name_len)    :: model_name

    model_name = "lfa"
    call CFG_add_get(cfg, "model%type", model_name, &
         "Which type of model is used")

    select case (model_name)
       case ("lfa")
          model_type = model_lfa
       case ("ee53")
          model_type = model_ee53
          model_has_energy_equation = .true.
    case default
       error stop "Unknown model (choices: lfa, ee53)"
    end select

  end subroutine model_initialize

end module m_model
