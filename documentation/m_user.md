# Customizing the code through `m_user.f90` {#m-user}

# Introduction

To provide more flexibility and to reduce the number of configuration parameters, certain methods can be customized by users. For example to use custom grid refinement or to write custom output.

## List of methods that can be modified

\snippet m_user_methods.f90 method_list

## How to set a custom method

In the `m_user.f90` file, first add your custom routine. For example, it could be a method for initial conditions called `my_init_cond`. This method can then be activated in the `user_initialize` subroutine:

    subroutine user_initialize(cfg, tree)
        type(CFG_t), intent(inout) :: cfg
        type(af_t), intent(inout) :: tree

        user_initial_conditions => my_init_cond
    end subroutine user_initialize

## Interfaces of the methods

For some of the interfaces, check out the [Afivo documentation](https://teunissen.net/afivo). The other interfaces are listed below.

\snippet m_user_methods.f90 interface_list

