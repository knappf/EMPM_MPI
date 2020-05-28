module int_arr

public :: n_nn,n_pp,n_pn

integer, dimension (:), allocatable,save  :: n_pp,n_nn,n_pn,n_pp_dd,n_nn_dd,n_pn_dd
double precision, dimension (:,:), allocatable :: v_pp,v_nn,v_pn,v_pp_dd,v_nn_dd,v_pn_dd
integer, dimension (:,:), allocatable :: icol_nn,icol_pp,icol_pn,icol_nn_dd,icol_pp_dd,icol_pn_dd
integer, dimension (:,:), allocatable :: irowc_nn,irowc_pp,irowc_pn,irowe_nn,irowe_pp,irowe_pn
integer, dimension (:,:), allocatable :: irowc_nn_dd,irowc_pp_dd,irowc_pn_dd,irowe_nn_dd,irowe_pp_dd,irowe_pn_dd

integer :: imax 


end module int_arr
