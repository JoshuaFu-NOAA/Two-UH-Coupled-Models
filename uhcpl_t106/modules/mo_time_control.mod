V26 mo_time_control
19 mo_time_control.f90 S582 0
10/14/2013  14:59:16
use mo_constants private
use mo_start_dataset private
use mo_control private
use mo_year private
use mo_constants private
use mo_start_dataset private
use mo_control private
use mo_year private
enduse
D 68 21 6 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
D 71 21 6 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
D 74 21 6 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
D 77 21 6 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
D 80 21 6 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
D 83 21 6 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
S 582 24 0 0 0 6 1 0 4658 10005 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 mo_time_control
S 584 23 0 0 0 9 715 582 4682 4 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 cd2dat
S 585 23 0 0 0 6 732 582 4689 4 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 im2day
S 587 23 0 0 0 9 749 582 4707 4 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 dtime
S 588 23 0 0 0 6 779 582 4713 4 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ntbase
S 589 23 0 0 0 6 778 582 4720 4 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ncbase
S 591 23 0 0 0 6 625 582 4744 4 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 nstep
S 592 23 0 0 0 6 663 582 4750 4 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ntimeadj
S 594 23 0 0 0 9 831 582 4772 4 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 dayl
S 596 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 625 6 2 mo_start_dataset nstep
R 663 6 40 mo_start_dataset ntimeadj
S 673 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 31 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 674 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 59 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 675 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 90 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 676 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 120 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 677 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 151 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 678 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 181 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 679 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 212 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 680 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 243 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 681 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 273 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 682 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 304 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 683 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 334 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 684 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 60 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 685 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 91 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 686 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 121 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 687 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 152 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 688 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 182 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 689 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 213 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 690 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 244 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 691 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 274 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 692 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 305 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 693 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 335 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 694 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 28 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 695 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 30 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 698 7 3 mo_year idays$ac
R 700 7 5 mo_year jdays$ac
R 702 7 7 mo_year kdays$ac
R 715 14 20 mo_year cd2dat
R 732 14 37 mo_year im2day
R 749 6 1 mo_control dtime
R 778 6 30 mo_control ncbase
R 779 6 31 mo_control ntbase
R 831 16 1 mo_constants dayl
S 867 16 0 0 0 6 1 582 6010 4 400000 A 0 0 0 0 0 0 0 0 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 time_inc_days
S 868 16 0 0 0 6 1 582 6024 4 400000 A 0 0 0 0 0 0 0 0 2 136 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 time_inc_hour
S 869 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 870 23 5 0 0 9 875 582 6038 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 time_control
S 871 1 3 1 0 6 1 870 6051 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 step_inc
S 872 1 3 1 0 6 1 870 6060 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 inc_type
S 873 1 3 1 0 6 1 870 6069 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 my_step
S 874 1 3 0 0 16 1 870 6077 4 1003000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 action
S 875 14 5 0 0 16 1 870 6038 4 1400000 A 0 0 0 0 0 0 0 44 3 0 0 874 0 0 0 0 0 0 0 0 0 36 0 582 0 0 0 0 time_control
F 875 3 871 872 873
S 876 23 5 0 0 9 880 582 6084 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 time_inc_sec
S 877 1 3 1 0 6 1 876 6051 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 step_inc
S 878 1 3 1 0 6 1 876 6060 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 inc_type
S 879 1 3 0 0 6 1 876 6097 4 1003000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 inc_sec
S 880 14 5 0 0 6 1 876 6084 4 1400000 A 0 0 0 0 0 0 0 48 2 0 0 879 0 0 0 0 0 0 0 0 0 96 0 582 0 0 0 0 time_inc_sec
F 880 2 877 878
S 881 23 5 0 0 9 885 582 6105 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 time_inc_steps
S 882 1 3 1 0 6 1 881 6051 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 step_inc
S 883 1 3 1 0 6 1 881 6060 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 inc_type
S 884 1 3 0 0 6 1 881 6120 4 1003000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 inc_steps
S 885 14 5 0 0 6 1 881 6105 4 1400000 A 0 0 0 0 0 0 0 51 2 0 0 884 0 0 0 0 0 0 0 0 0 125 0 582 0 0 0 0 time_inc_steps
F 885 2 882 883
A 12 2 0 0 0 6 596 0 0 0 12 0 0 0 0 0 0 0 0 0
A 14 2 0 0 0 6 673 0 0 0 14 0 0 0 0 0 0 0 0 0
A 15 2 0 0 0 6 674 0 0 0 15 0 0 0 0 0 0 0 0 0
A 16 2 0 0 0 6 675 0 0 0 16 0 0 0 0 0 0 0 0 0
A 17 2 0 0 0 6 676 0 0 0 17 0 0 0 0 0 0 0 0 0
A 18 2 0 0 0 6 677 0 0 0 18 0 0 0 0 0 0 0 0 0
A 19 2 0 0 0 6 678 0 0 0 19 0 0 0 0 0 0 0 0 0
A 20 2 0 0 0 6 679 0 0 0 20 0 0 0 0 0 0 0 0 0
A 21 2 0 0 0 6 680 0 0 0 21 0 0 0 0 0 0 0 0 0
A 22 2 0 0 0 6 681 0 0 0 22 0 0 0 0 0 0 0 0 0
A 23 2 0 0 0 6 682 0 0 0 23 0 0 0 0 0 0 0 0 0
A 24 2 0 0 0 6 683 0 0 0 24 0 0 0 0 0 0 0 0 0
A 38 2 0 0 0 6 684 0 0 0 38 0 0 0 0 0 0 0 0 0
A 39 2 0 0 0 6 685 0 0 0 39 0 0 0 0 0 0 0 0 0
A 40 2 0 0 0 6 686 0 0 0 40 0 0 0 0 0 0 0 0 0
A 41 2 0 0 0 6 687 0 0 0 41 0 0 0 0 0 0 0 0 0
A 42 2 0 0 0 6 688 0 0 0 42 0 0 0 0 0 0 0 0 0
A 43 2 0 0 0 6 689 0 0 0 43 0 0 0 0 0 0 0 0 0
A 44 2 0 0 0 6 690 0 0 0 44 0 0 0 0 0 0 0 0 0
A 45 2 0 0 0 6 691 0 0 0 45 0 0 0 0 0 0 0 0 0
A 46 2 0 0 0 6 692 0 0 0 46 0 0 0 0 0 0 0 0 0
A 47 2 0 0 0 6 693 0 0 0 47 0 0 0 0 0 0 0 0 0
A 61 2 0 0 0 6 694 0 0 0 61 0 0 0 0 0 0 0 0 0
A 62 2 0 0 0 6 695 0 0 0 62 0 0 0 0 0 0 0 0 0
A 89 1 0 1 0 68 698 0 0 0 0 0 0 0 0 0 0 0 0 0
A 103 1 0 1 0 74 700 0 0 0 0 0 0 0 0 0 0 0 0 0
A 117 1 0 1 0 80 702 0 0 0 0 0 0 0 0 0 0 0 0 0
A 136 2 0 0 0 6 869 0 0 0 136 0 0 0 0 0 0 0 0 0
Z
J 10 1 1
V 89 68 7 0
R 0 71 0 0
A 0 6 0 0 1 2 1
A 0 6 0 0 1 14 1
A 0 6 0 0 1 15 1
A 0 6 0 0 1 16 1
A 0 6 0 0 1 17 1
A 0 6 0 0 1 18 1
A 0 6 0 0 1 19 1
A 0 6 0 0 1 20 1
A 0 6 0 0 1 21 1
A 0 6 0 0 1 22 1
A 0 6 0 0 1 23 1
A 0 6 0 0 1 24 0
J 11 1 1
V 103 74 7 0
R 0 77 0 0
A 0 6 0 0 1 2 1
A 0 6 0 0 1 14 1
A 0 6 0 0 1 38 1
A 0 6 0 0 1 39 1
A 0 6 0 0 1 40 1
A 0 6 0 0 1 41 1
A 0 6 0 0 1 42 1
A 0 6 0 0 1 43 1
A 0 6 0 0 1 44 1
A 0 6 0 0 1 45 1
A 0 6 0 0 1 46 1
A 0 6 0 0 1 47 0
J 12 1 1
V 117 80 7 0
R 0 83 0 0
A 0 6 0 0 1 14 1
A 0 6 0 0 1 61 1
A 0 6 0 0 1 14 1
A 0 6 0 0 1 62 1
A 0 6 0 0 1 14 1
A 0 6 0 0 1 62 1
A 0 6 0 0 1 14 1
A 0 6 0 0 1 14 1
A 0 6 0 0 1 62 1
A 0 6 0 0 1 14 1
A 0 6 0 0 1 62 1
A 0 6 0 0 1 14 0
Z