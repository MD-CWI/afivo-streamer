!> This module contains routines for restriction: going from fine to coarse
!> variables.
module m_a$D_restrict

  use m_a$D_types

  implicit none
  private

  public :: a$D_restrict_to_box
  public :: a$D_restrict_to_boxes
  public :: a$D_restrict_tree
  public :: a$D_restrict_box
  public :: a$D_restrict_box_cont
  public :: a$D_restrict_box_fc

contains

  !> Restrict the children of a box to the box (e.g., in 2D, average the values
  !> at the four children to get the value for the parent)
  subroutine a$D_restrict_to_box(boxes, id, iv, iv_to, i_eps, med, s_flag)
    type(box$D_t), intent(inout)   :: boxes(:) !< List of all the boxes
    integer, intent(in)            :: id       !< Box whose children will be restricted to it
    integer, intent(in)            :: iv       !< Variable to restrict
    integer, intent(in), optional  :: iv_to    !< Destination (if /= iv)
    integer, intent(in), optional  :: i_eps
    real(dp), intent(in), optional :: med
    logical, intent(in), optional  :: s_flag
    integer                        :: nc, i_c, c_id
    logical                        :: surf
    
    surf = .false.
    if(present(s_flag)) surf = s_flag

    nc = boxes(id)%n_cell
    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == af_no_box) cycle
       if (surf) then
         call a$D_restrict_box_fc(boxes(c_id), boxes(id), iv, iv_to, i_eps, med)
       else
         call a$D_restrict_box(boxes(c_id), boxes(id), iv, iv_to, i_eps, med)
       end if
    end do
  end subroutine a$D_restrict_to_box

  !> Restrict the children of boxes ids(:) to them.
  subroutine a$D_restrict_to_boxes(boxes, ids, iv, iv_to, i_eps, med, s_flag)
    type(box$D_t), intent(inout)   :: boxes(:) !< List of all the boxes
    integer, intent(in)            :: ids(:)   !< Boxes whose children will be restricted to it
    integer, intent(in)            :: iv       !< Variable to restrict
    integer, intent(in), optional  :: iv_to    !< Destination (if /= iv)
    logical, intent(in), optional  :: s_flag
    integer, intent(in), optional  :: i_eps
    real(dp), intent(in), optional :: med
    integer                        :: i

    !$omp parallel do
    do i = 1, size(ids)
       call a$D_restrict_to_box(boxes, ids(i), iv, iv_to, i_eps, med, s_flag)
    end do
    !$omp end parallel do
  end subroutine a$D_restrict_to_boxes

  !> Restrict variables iv to all parent boxes, from the highest to the lowest level
  subroutine a$D_restrict_tree(tree, iv, iv_to, i_eps, med, s_flag)
    type(a$D_t), intent(inout)     :: tree  !< Tree to restrict on
    integer, intent(in)            :: iv    !< Variable to restrict
    integer, intent(in), optional  :: iv_to !< Destination (if /= iv)
    logical, intent(in), optional  :: s_flag
    integer, intent(in), optional  :: i_eps
    real(dp), intent(in), optional :: med
    integer                        :: lvl

    if (.not. tree%ready) stop "Tree not ready"
    do lvl = tree%highest_lvl-1, lbound(tree%lvls, 1), -1
       call a$D_restrict_to_boxes(tree%boxes, tree%lvls(lvl)%parents, iv, iv_to, i_eps, med, s_flag)
    end do
  end subroutine a$D_restrict_tree

  !> Restriction of child box (box_c) to its parent (box_p)
  subroutine a$D_restrict_box(box_c, box_p, iv, iv_to, i_eps, med)
    use m_a$D_utils, only:a$D_heaviside 
    type(box$D_t), intent(in)       :: box_c         !< Child box to restrict
    type(box$D_t), intent(inout)    :: box_p         !< Parent box to restrict to
    integer, intent(in)             :: iv            !< Variable to restrict
    integer, intent(in), optional   :: iv_to         !< Destination (if /= iv) 
    integer, intent(in), optional   :: i_eps
    real(dp), intent(in), optional  :: med
    integer                       :: i, j, i_f, j_f, i_c, j_c, i_dest
    integer                       :: hnc, ix_offset($D)
#if $D == 2
    real(dp)                      :: w1, w2
    real(dp)                      :: in_out(0:1, 0:1)
#elif $D == 3
    integer                       :: k, k_f, k_c
#endif

    hnc       = ishft(box_c%n_cell, -1) ! n_cell / 2
    ix_offset = a$D_get_child_offset(box_c)

    if (present(iv_to)) then
       i_dest = iv_to
    else
       i_dest = iv
    end if

#if $D == 2
    if (box_p%coord_t == af_cyl) then
       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1
             
             in_out(0:1, 0:1) = a$D_heaviside(box_c%cc(i_f:i_f+1, j_f:j_f+1, i_eps), med) 

             call a2_cyl_child_weights(box_p, i, w1, w2)
             box_p%cc(i_c, j_c, i_dest) = (w1 * sum(in_out(0, 0:1)*box_c%cc(i_f, j_f:j_f+1, iv)) + &
                  w2 * sum(in_out(1, 0:1)*box_c%cc(i_f+1, j_f:j_f+1, iv))) / sum(in_out(:,:))
          end do
       end do
    else
       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1
             
             in_out(0:1, 0:1) = a$D_heaviside(box_c%cc(i_f:i_f+1, j_f:j_f+1, i_eps), med)
             if (sum(in_out(:,:)) > epsilon(1.0_dp)) then
               if (box_p%cc(i_c, j_c, i_eps) <= med) then
                 box_p%cc(i_c, j_c, i_dest) = sum(in_out(0:1, 0:1)*box_c%cc(i_f:i_f+1, j_f:j_f+1, iv)) / sum(in_out(:,:))
               else
                 box_p%cc(i_c, j_c, i_dest) = sum((1-in_out(0:1, 0:1))*box_c%cc(i_f:i_f+1, j_f:j_f+1, iv)) / sum(1-in_out(:,:))
               end if
             else
               box_p%cc(i_c, j_c, i_dest) = 0.25_dp*sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, iv))
             end if

          end do
       end do
    endif
#elif $D == 3
    do k = 1, hnc
       k_c = ix_offset(3) + k
       k_f = 2 * k - 1
       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1
             box_p%cc(i_c, j_c, k_c, i_dest) = 0.125_dp * &
                  sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, k_f:k_f+1, iv))
          end do
       end do
    end do
#endif
  end subroutine a$D_restrict_box
 
 
  subroutine a$D_restrict_box_cont(box_c, box_p, iv, iv_to, i_eps)
    use m_a$D_utils, only:a$D_heaviside 
    type(box$D_t), intent(in)       :: box_c         !< Child box to restrict
    type(box$D_t), intent(inout)    :: box_p         !< Parent box to restrict to
    integer, intent(in)             :: iv            !< Variable to restrict
    integer, intent(in), optional   :: iv_to         !< Destination (if /= iv) 
    integer, intent(in), optional   :: i_eps
    integer                       :: i, j, i_f, j_f, i_c, j_c, i_dest
    integer                       :: hnc, ix_offset($D)
#if $D == 2
    real(dp)                      :: w1, w2
    real(dp)                      :: inv_eps(0:1, 0:1)
#elif $D == 3
    integer                       :: k, k_f, k_c
#endif

    hnc       = ishft(box_c%n_cell, -1) ! n_cell / 2
    ix_offset = a$D_get_child_offset(box_c)

    if (present(iv_to)) then
       i_dest = iv_to
    else
       i_dest = iv
    end if

#if $D == 2
    if (box_p%coord_t == af_cyl) then
       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1
             
             inv_eps(0:1, 0:1) = 1.0_dp/box_c%cc(i_f:i_f+1, j_f:j_f+1, i_eps)

             call a2_cyl_child_weights(box_p, i, w1, w2)
             box_p%cc(i_c, j_c, i_dest) = (w1 * sum(inv_eps(0, 0:1)*box_c%cc(i_f, j_f:j_f+1, iv)) + &
                  w2 * sum(inv_eps(1, 0:1)*box_c%cc(i_f+1, j_f:j_f+1, iv))) / sum(inv_eps(:,:))
          end do
       end do
    else
       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1
             inv_eps(0:1, 0:1) = (box_p%cc(i_c, j_c, i_eps)+box_c%cc(i_f:i_f+1, j_f:j_f+1, i_eps))/&
                                 (box_p%cc(i_c, j_c, i_eps)*box_c%cc(i_f:i_f+1, j_f:j_f+1, i_eps))
             
             box_p%cc(i_c, j_c, i_dest) = sum(inv_eps(0:1, 0:1)*box_c%cc(i_f:i_f+1, j_f:j_f+1, iv)) / sum(inv_eps(:,:))
          end do
       end do
    endif
#elif $D == 3
    do k = 1, hnc
       k_c = ix_offset(3) + k
       k_f = 2 * k - 1
       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1
             box_p%cc(i_c, j_c, k_c, i_dest) = 0.125_dp * &
                  sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, k_f:k_f+1, iv))
          end do
       end do
    end do
#endif
  end subroutine a$D_restrict_box_cont
  
  
  !> Restriction of child box (box_c) to its parent (box_p)
  subroutine a$D_restrict_box_fc(box_c, box_p, s_iv, s_iv_to, i_eps, med)
    use m_a$D_utils, only:a$D_border
    type(box$D_t), intent(in)       :: box_c         !< Child box to restrict
    type(box$D_t), intent(inout)    :: box_p         !< Parent box to restrict to
    integer, intent(in)             :: s_iv            !< Variable to restrict
    integer, intent(in), optional   :: s_iv_to         !< Destination (if /= iv) 
    integer, intent(in), optional   :: i_eps
    real(dp), intent(in), optional  :: med
    integer                       :: i, j, i_f, j_f, i_c, j_c, i_dest
    integer                       :: hnc, ix_offset($D), ii
#if $D == 2
    real(dp)                      :: w1, w2
#elif $D == 3
    integer                       :: k, k_f, k_c
#endif

    hnc       = ishft(box_c%n_cell, -1) ! n_cell / 2
    ix_offset = a$D_get_child_offset(box_c)

    if (present(s_iv_to)) then
       i_dest = s_iv_to
    else
       i_dest = s_iv
    end if

#if $D == 2
    if (box_p%coord_t == af_cyl) then
       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1
             
          end do
       end do
    else
       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1
              
             if (a$D_border(box_p%cc(i_c, j_c, i_eps), box_p%cc(i_c-1, j_c, i_eps), med)) then    
               if (i_f == 1) then
                 box_p%fc(i_c, j_c, 1, i_dest) = 0.5_dp*sum(box_c%fc(i_f:i_f+1, j_f, 1, s_iv)) + &
                                               0.5_dp*sum(box_c%fc(i_f:i_f+1, j_f+1, 1, s_iv))
               else       
                 box_p%fc(i_c, j_c, 1, i_dest) = 0.5_dp*sum(box_c%fc(i_f-1:i_f+1, j_f, 1, s_iv)) + &
                                               0.5_dp*sum(box_c%fc(i_f-1:i_f+1, j_f+1, 1, s_iv))
               end if
             end if
               
             if (a$D_border(box_p%cc(i_c, j_c, i_eps), box_p%cc(i_c, j_c-1, i_eps), med)) then    
               if (j_f == 1) then
                 box_p%fc(i_c, j_c, 2, i_dest) = 0.5_dp*sum(box_c%fc(i_f, j_f:j_f+1, 2, s_iv)) + &
                                               0.5_dp*sum(box_c%fc(i_f+1, j_f:j_f+1, 2, s_iv))
               else       
                 box_p%fc(i_c, j_c, 2, i_dest) = 0.5_dp*sum(box_c%fc(i_f, j_f-1:j_f+1, 2, s_iv)) + &
                                               0.5_dp*sum(box_c%fc(i_f+1, j_f-1:j_f+1, 2, s_iv))
               end if
             end if
             
          end do
       end do
    endif
#elif $D == 3
    do k = 1, hnc
       k_c = ix_offset(3) + k
       k_f = 2 * k - 1
       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1

          end do
       end do
    end do
#endif
  end subroutine a$D_restrict_box_fc

end module m_a$D_restrict
