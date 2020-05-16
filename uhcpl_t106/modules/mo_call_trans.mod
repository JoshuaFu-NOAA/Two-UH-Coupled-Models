V26 mo_call_trans
17 mo_call_trans.f90 S582 0
10/14/2013  14:59:18
use mo_linked_list private
use mo_test_trans private
use mo_doctor private
use mo_mpi private
use mo_decomposition private
use mo_linked_list private
use mo_test_trans private
use mo_doctor private
use mo_mpi private
use mo_decomposition private
enduse
D 1059 24 4090 272 4080 3
D 1084 18 23
D 1086 18 212
D 1088 21 6 1 3 189 0 0 0 0 0
 0 189 3 3 189 189
D 1144 24 4126 24 4082 7
S 582 24 0 0 0 6 1 0 4658 10015 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 mo_call_trans
S 584 23 0 0 0 9 649 582 4689 10 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 dc
S 586 23 0 0 0 9 650 582 4713 10 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 dcl
S 588 23 0 0 0 9 651 582 4737 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 debug_parallel
S 589 23 0 0 0 9 653 582 4752 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 any_col_1d
S 591 23 0 0 0 9 1059 582 4770 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 p_io
S 593 23 0 0 0 6 608 582 4785 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 nerr
S 595 19 0 0 0 9 1 582 4804 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 372 2 0 0 0 0 0 582 0 0 0 0 test_spectral
O 595 2 3593 3607
S 596 19 0 0 0 9 1 582 4818 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 369 3 0 0 0 0 0 582 0 0 0 0 test_legendre
O 596 3 3618 3635 3649
S 597 19 0 0 0 9 1 582 4832 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 365 3 0 0 0 0 0 582 0 0 0 0 test_gridpoint
O 597 3 3688 3705 3719
S 598 19 0 0 0 9 1 582 4847 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 361 2 0 0 0 0 0 582 0 0 0 0 test_symasym
O 598 2 3660 3677
S 599 23 0 0 0 9 3731 582 4860 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 test_zonmean
S 600 19 0 0 0 9 1 582 4873 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 357 4 0 0 0 0 0 582 0 0 0 0 test_row
O 600 4 3771 3758 3748 3787
R 608 16 3 mo_doctor nerr
S 641 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 649 7 3 mo_decomposition global_decomposition
R 650 6 4 mo_decomposition local_decomposition
R 651 6 5 mo_decomposition debug_parallel
R 653 6 7 mo_decomposition any_col_1d
S 877 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 887 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 956 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 128 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 1059 6 19 mo_mpi p_io
R 3593 14 48 mo_test_trans test_spectral3
R 3607 14 62 mo_test_trans test_spectral2
R 3618 14 73 mo_test_trans test_legendre4
R 3635 14 90 mo_test_trans test_legendre3
R 3649 14 104 mo_test_trans test_legendre0
R 3660 14 115 mo_test_trans test_symasym4
R 3677 14 132 mo_test_trans test_symasym2
R 3688 14 143 mo_test_trans test_gridpoint4
R 3705 14 160 mo_test_trans test_gridpoint3
R 3719 14 174 mo_test_trans test_gridpoint2
R 3731 14 186 mo_test_trans test_zonmean
R 3748 14 203 mo_test_trans test_row1
R 3758 14 213 mo_test_trans test_row2
R 3771 14 226 mo_test_trans test_row3
R 3787 14 242 mo_test_trans test_row0
S 4035 27 0 0 0 9 4047 582 19321 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 spectral_to_legendre
S 4036 27 0 0 0 6 4049 582 19342 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 legendre_to_fourier
S 4037 27 0 0 0 9 4053 582 19362 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 fourier_to_gridpoint
S 4038 27 0 0 0 9 4055 582 19383 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 gridpoint_to_fourier
S 4039 27 0 0 0 9 4051 582 19404 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 fourier_to_legendre
S 4040 27 0 0 0 6 4045 582 19424 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 legendre_to_spectral
S 4041 27 0 0 0 9 4060 582 19445 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 test_memory_f
S 4042 27 0 0 0 9 4176 582 19459 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 test_memory_gp
S 4043 27 0 0 0 9 4057 582 19474 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 test_scan_buffer
S 4044 27 0 0 0 9 4180 582 19491 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 test_row_buffer
S 4045 23 5 0 0 0 4046 582 19424 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 legendre_to_spectral
S 4046 14 5 0 0 0 1 4045 19424 0 400000 A 0 0 0 0 0 0 0 763 0 0 0 0 0 0 0 0 0 0 0 0 0 53 0 582 0 0 0 0 legendre_to_spectral
F 4046 0
S 4047 23 5 0 0 0 4048 582 19321 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 spectral_to_legendre
S 4048 14 5 0 0 0 1 4047 19321 0 400000 A 0 0 0 0 0 0 0 764 0 0 0 0 0 0 0 0 0 0 0 0 0 85 0 582 0 0 0 0 spectral_to_legendre
F 4048 0
S 4049 23 5 0 0 0 4050 582 19342 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 legendre_to_fourier
S 4050 14 5 0 0 0 1 4049 19342 0 400000 A 0 0 0 0 0 0 0 765 0 0 0 0 0 0 0 0 0 0 0 0 0 117 0 582 0 0 0 0 legendre_to_fourier
F 4050 0
S 4051 23 5 0 0 0 4052 582 19404 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 fourier_to_legendre
S 4052 14 5 0 0 0 1 4051 19404 0 400000 A 0 0 0 0 0 0 0 766 0 0 0 0 0 0 0 0 0 0 0 0 0 135 0 582 0 0 0 0 fourier_to_legendre
F 4052 0
S 4053 23 5 0 0 0 4054 582 19362 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 fourier_to_gridpoint
S 4054 14 5 0 0 0 1 4053 19362 0 400000 A 0 0 0 0 0 0 0 767 0 0 0 0 0 0 0 0 0 0 0 0 0 154 0 582 0 0 0 0 fourier_to_gridpoint
F 4054 0
S 4055 23 5 0 0 0 4056 582 19383 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gridpoint_to_fourier
S 4056 14 5 0 0 0 1 4055 19383 0 400000 A 0 0 0 0 0 0 0 768 0 0 0 0 0 0 0 0 0 0 0 0 0 197 0 582 0 0 0 0 gridpoint_to_fourier
F 4056 0
S 4057 23 5 0 0 0 4059 582 19474 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 test_scan_buffer
S 4058 1 3 1 0 28 1 4057 5071 14 43000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 text
S 4059 14 5 0 0 0 1 4057 19474 0 400000 A 0 0 0 0 0 0 0 769 1 0 0 0 0 0 0 0 0 0 0 0 0 215 0 582 0 0 0 0 test_scan_buffer
F 4059 1 4058
S 4060 23 5 0 0 0 4062 582 19445 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 test_memory_f
S 4061 1 3 1 0 28 1 4060 5071 14 43000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 text
S 4062 14 5 0 0 0 1 4060 19445 0 400000 A 0 0 0 0 0 0 0 771 1 0 0 0 0 0 0 0 0 0 0 0 0 254 0 582 0 0 0 0 test_memory_f
F 4062 1 4061
S 4070 3 0 0 0 1084 1 1 0 0 0 A 0 0 0 0 0 0 0 0 19542 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 8 20 20 20 20 20 20 20 20
S 4071 3 0 0 0 1086 1 1 0 0 0 A 0 0 0 0 0 0 0 0 19551 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 9 75 6e 64 65 66 69 6e 65 64
S 4072 3 0 0 0 16 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16
R 4080 25 5 mo_linked_list memory_info
R 4082 25 7 mo_linked_list list
R 4090 5 15 mo_linked_list name memory_info
R 4091 5 16 mo_linked_list dim_1 memory_info
R 4092 5 17 mo_linked_list dim_2 memory_info
R 4093 5 18 mo_linked_list dim_3 memory_info
R 4094 5 19 mo_linked_list dim_4 memory_info
R 4095 5 20 mo_linked_list gdim_1 memory_info
R 4096 5 21 mo_linked_list gdim_2 memory_info
R 4097 5 22 mo_linked_list gdim_3 memory_info
R 4098 5 23 mo_linked_list gdim_4 memory_info
R 4099 5 24 mo_linked_list ndim memory_info
R 4100 5 25 mo_linked_list gribtable memory_info
R 4101 5 26 mo_linked_list gribcode memory_info
R 4102 5 27 mo_linked_list io_var_indx memory_info
R 4103 5 28 mo_linked_list io_var_id memory_info
R 4104 5 29 mo_linked_list io_name memory_info
R 4105 5 30 mo_linked_list io_unit memory_info
R 4106 5 31 mo_linked_list outint memory_info
R 4107 5 32 mo_linked_list accumulate memory_info
R 4108 5 33 mo_linked_list assign memory_info
R 4109 5 34 mo_linked_list restart memory_info
R 4120 6 45 mo_linked_list empty_info$ac
R 4126 5 51 mo_linked_list first_list_element list
R 4128 5 53 mo_linked_list first_list_element$p list
R 4130 5 55 mo_linked_list memory_used list
R 4131 5 56 mo_linked_list list_elements list
S 4176 23 5 0 0 0 4179 582 19459 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 test_memory_gp
S 4177 1 3 1 0 1144 1 4176 17633 14 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gp
S 4178 1 3 1 0 28 1 4176 5071 14 43000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 text
S 4179 14 5 0 0 0 1 4176 19459 0 400000 A 0 0 0 0 0 0 0 796 2 0 0 0 0 0 0 0 0 0 0 0 0 273 0 582 0 0 0 0 test_memory_gp
F 4179 2 4177 4178
S 4180 23 5 0 0 0 4183 582 19491 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 test_row_buffer
S 4181 1 3 1 0 6 1 4180 17750 14 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j
S 4182 1 3 1 0 28 1 4180 5071 14 43000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 text
S 4183 14 5 0 0 0 1 4180 19491 0 400000 A 0 0 0 0 0 0 0 799 2 0 0 0 0 0 0 0 0 0 0 0 0 292 0 582 0 0 0 0 test_row_buffer
F 4183 2 4181 4182
A 23 2 0 0 0 6 641 0 0 0 23 0 0 0 0 0 0 0 0 0
A 189 2 0 0 0 6 877 0 0 0 189 0 0 0 0 0 0 0 0 0
A 212 2 0 0 0 6 887 0 0 0 212 0 0 0 0 0 0 0 0 0
A 563 2 0 0 190 6 956 0 0 0 563 0 0 0 0 0 0 0 0 0
A 2939 2 0 0 2704 16 4072 0 0 0 2939 0 0 0 0 0 0 0 0 0
A 2984 2 0 0 2627 1084 4070 0 0 0 2984 0 0 0 0 0 0 0 0 0
A 2985 2 0 0 2637 1086 4071 0 0 0 2985 0 0 0 0 0 0 0 0 0
A 3011 1 0 0 2286 1059 4120 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 60 1 1
V 3011 1059 7 0
S 0 1059 0 0 0
A 0 1084 0 0 1 2984 1
A 0 6 0 0 1 2 1
A 0 6 0 0 1 3 1
A 0 6 0 0 1 3 1
A 0 6 0 0 1 3 1
A 0 6 0 0 1 2 1
A 0 6 0 0 1 2 1
A 0 6 0 0 1 2 1
A 0 6 0 0 1 2 1
A 0 6 0 0 1 2 1
A 0 6 0 0 1 563 1
A 0 6 0 0 1 2 1
R 0 1088 0 1
A 0 6 0 0 1 2 1
A 0 6 0 0 1 2 1
A 0 6 0 0 1 2 1
A 0 6 0 0 1 2 0
A 0 6 0 0 1 2 1
A 0 1086 0 0 1 2985 1
A 0 1086 0 0 1 2985 1
A 0 6 0 0 1 2 1
A 0 16 0 0 1 2939 1
A 0 16 0 0 1 2939 1
A 0 16 0 0 1 2939 0
Z
