V26 mo_truncation
17 mo_truncation.f90 S582 0
10/14/2013  14:59:15
use mo_doctor private
use mo_parameters private
use mo_doctor private
use mo_parameters private
enduse
D 58 21 6 1 55 53 0 1 0 0 1
 45 49 51 45 49 47
D 61 21 6 1 0 42 0 0 0 0 0
 0 42 0 3 42 0
D 64 21 6 1 64 63 0 1 0 0 1
 58 61 62 58 61 59
D 67 21 6 1 0 42 0 0 0 0 0
 0 42 0 3 42 0
D 70 21 6 1 3 21 0 0 0 0 0
 0 21 3 3 21 21
D 73 21 6 1 3 29 0 0 0 0 0
 0 29 3 3 29 29
D 76 21 6 1 3 29 0 0 0 0 0
 0 29 3 3 29 29
D 79 21 9 1 3 29 0 0 0 0 0
 0 29 3 3 29 29
S 582 24 0 0 0 6 1 0 4658 10015 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 mo_truncation
S 584 23 0 0 0 6 609 582 4686 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 jphgl
S 585 23 0 0 0 6 608 582 4692 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 jpmp1
S 586 23 0 0 0 6 604 582 4698 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 jpnlev
S 588 23 0 0 0 6 615 582 4715 14 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 nout
S 592 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 999 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 596 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 107 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 604 16 6 mo_parameters jpnlev
R 608 16 10 mo_parameters jpmp1
R 609 16 11 mo_parameters jphgl
R 615 16 2 mo_doctor nout
S 629 7 4 0 4 73 630 582 4881 800004 100 A 0 0 0 0 0 0 0 0 0 0 0 0 652 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 nmp
S 630 7 4 0 4 76 632 582 4885 800004 100 A 0 0 0 0 0 432 0 0 0 0 0 0 652 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 nnp
S 631 7 4 0 4 79 1 582 4889 800004 100 A 0 0 0 0 0 0 0 0 0 0 0 0 653 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 am
S 632 7 4 0 4 70 1 582 4892 800004 100 A 0 0 0 0 0 864 0 0 0 0 0 0 652 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ntrn
S 633 7 6 0 0 58 1 582 4897 10a00004 51 A 0 0 0 0 0 0 639 0 0 0 641 0 0 0 0 0 0 0 0 638 0 0 640 582 0 0 0 0 ntrm
S 634 7 6 0 0 64 1 582 4902 10a00004 51 A 0 0 0 0 0 0 649 0 0 0 651 0 0 0 0 0 0 0 0 648 0 0 650 582 0 0 0 0 ntrk
S 635 27 0 0 0 9 655 582 4907 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 scpar
S 636 6 4 0 0 6 647 582 4913 40800016 0 A 0 0 0 0 0 0 0 0 0 0 0 0 654 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 z_b_0
S 637 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 0 0 0 0 0 18 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 638 8 4 0 0 61 649 582 4919 40822004 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 652 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ntrm$sd
S 639 6 4 0 0 7 640 582 4927 40802001 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 652 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ntrm$p
S 640 6 4 0 0 7 638 582 4934 40802000 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 652 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ntrm$o
S 641 22 1 0 0 6 1 582 4941 40000000 1000 A 0 0 0 0 0 0 0 633 0 0 0 0 638 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ntrm$arrdsc
S 642 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 0 0 0 0 0 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 643 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 0 0 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 644 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 0 0 0 0 0 17 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 645 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 646 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 0 0 0 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 647 6 4 0 0 6 1 582 4953 40800016 0 A 0 0 0 0 0 4 0 0 0 0 0 0 654 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 z_b_1
S 648 8 4 0 0 67 629 582 4959 40822004 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 652 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ntrk$sd
S 649 6 4 0 0 7 650 582 4967 40802001 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 652 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ntrk$p
S 650 6 4 0 0 7 648 582 4974 40802000 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 652 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ntrk$o
S 651 22 1 0 0 6 1 582 4981 40000000 1000 A 0 0 0 0 0 0 0 634 0 0 0 0 648 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ntrk$arrdsc
S 652 11 0 0 4 9 628 582 4993 40800000 801000 A 0 0 0 0 0 5036 0 0 639 632 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _mo_truncation$0
S 653 11 0 0 4 9 652 582 5010 40800000 801000 A 0 0 0 0 0 856 0 0 631 631 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _mo_truncation$2
S 654 11 0 0 0 9 653 582 5027 40800010 801000 A 0 0 0 0 0 8 0 0 636 647 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _mo_truncation$4
S 655 23 5 0 0 0 659 582 4907 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 scpar
S 656 1 3 1 0 6 1 655 5044 14 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nm
S 657 1 3 1 0 6 1 655 5047 14 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nn
S 658 1 3 1 0 6 1 655 5050 14 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nk
S 659 14 5 0 0 0 1 655 4907 0 400000 A 0 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 0 0 28 0 582 0 0 0 0 scpar
F 659 3 656 657 658
A 21 2 0 0 0 6 592 0 0 0 21 0 0 0 0 0 0 0 0 0
A 29 2 0 0 0 6 596 0 0 0 29 0 0 0 0 0 0 0 0 0
A 42 2 0 0 0 6 637 0 0 0 42 0 0 0 0 0 0 0 0 0
A 43 2 0 0 0 6 642 0 0 0 43 0 0 0 0 0 0 0 0 0
A 44 1 0 1 0 61 638 0 0 0 0 0 0 0 0 0 0 0 0 0
A 45 10 0 0 0 6 44 1 0 0 0 0 0 0 0 0 0 0 0 0
X 1 43
A 46 2 0 0 0 6 643 0 0 0 46 0 0 0 0 0 0 0 0 0
A 47 10 0 0 45 6 44 4 0 0 0 0 0 0 0 0 0 0 0 0
X 1 46
A 48 4 0 0 0 6 47 0 3 0 0 0 0 2 0 0 0 0 0 0
A 49 4 0 0 0 6 45 0 48 0 0 0 0 1 0 0 0 0 0 0
A 50 2 0 0 0 6 644 0 0 0 50 0 0 0 0 0 0 0 0 0
A 51 10 0 0 47 6 44 7 0 0 0 0 0 0 0 0 0 0 0 0
X 1 50
A 52 2 0 0 0 6 645 0 0 0 52 0 0 0 0 0 0 0 0 0
A 53 10 0 0 51 6 44 10 0 0 0 0 0 0 0 0 0 0 0 0
X 1 52
A 54 2 0 0 0 6 646 0 0 0 54 0 0 0 0 0 0 0 0 0
A 55 10 0 0 53 6 44 13 0 0 0 0 0 0 0 0 0 0 0 0
X 1 54
A 57 1 0 1 0 67 648 0 0 0 0 0 0 0 0 0 0 0 0 0
A 58 10 0 0 0 6 57 1 0 0 0 0 0 0 0 0 0 0 0 0
X 1 43
A 59 10 0 0 58 6 57 4 0 0 0 0 0 0 0 0 0 0 0 0
X 1 46
A 60 4 0 0 0 6 59 0 3 0 0 0 0 2 0 0 0 0 0 0
A 61 4 0 0 0 6 58 0 60 0 0 0 0 1 0 0 0 0 0 0
A 62 10 0 0 59 6 57 7 0 0 0 0 0 0 0 0 0 0 0 0
X 1 50
A 63 10 0 0 62 6 57 10 0 0 0 0 0 0 0 0 0 0 0 0
X 1 52
A 64 10 0 0 63 6 57 13 0 0 0 0 0 0 0 0 0 0 0 0
X 1 54
Z
Z