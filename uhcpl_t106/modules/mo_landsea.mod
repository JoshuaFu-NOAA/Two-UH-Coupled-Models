V26 mo_landsea
14 mo_landsea.f90 S582 0
10/14/2013  14:59:14
use mo_doctor private
use mo_control private
use mo_doctor private
use mo_control private
enduse
D 73 21 6 2 54 53 0 1 0 0 1
 40 43 50 40 43 41
 45 49 52 45 49 47
D 76 21 6 1 0 38 0 0 0 0 0
 0 38 0 3 38 0
D 79 21 6 2 58 57 0 0 0 0 0
 0 55 3 3 55 55
 0 56 55 3 56 56
D 82 21 6 2 73 72 0 1 0 0 1
 62 65 70 62 65 63
 66 69 71 66 69 67
D 85 21 6 1 0 38 0 0 0 0 0
 0 38 0 3 38 0
D 88 21 6 2 88 87 0 1 0 0 1
 77 80 85 77 80 78
 81 84 86 81 84 82
D 91 21 6 1 0 38 0 0 0 0 0
 0 38 0 3 38 0
D 94 21 9 1 97 96 0 1 0 0 1
 91 94 95 91 94 92
D 97 21 6 1 0 12 0 0 0 0 0
 0 12 0 3 12 0
D 100 21 9 1 106 105 0 1 0 0 1
 100 103 104 100 103 101
D 103 21 6 1 0 12 0 0 0 0 0
 0 12 0 3 12 0
D 106 21 9 1 3 55 0 0 0 0 0
 0 55 3 3 55 55
D 109 21 9 1 3 56 0 0 0 0 0
 0 56 3 3 56 56
D 112 21 9 2 121 120 0 1 0 0 1
 110 113 118 110 113 111
 114 117 119 114 117 115
D 115 21 6 1 0 38 0 0 0 0 0
 0 38 0 3 38 0
D 118 21 9 1 3 123 0 0 1 0 0
 0 122 3 3 123 123
D 121 21 9 1 3 125 0 0 1 0 0
 0 124 3 3 125 125
D 124 21 9 2 126 132 1 1 0 0 1
 3 127 3 3 127 128
 3 129 130 3 129 131
D 127 21 6 2 133 138 0 0 1 0 0
 0 134 3 3 135 135
 0 136 135 3 137 137
D 130 21 6 2 58 57 0 0 0 0 0
 0 55 3 3 55 55
 0 56 55 3 56 56
S 582 24 0 0 0 6 1 0 4658 10005 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 mo_landsea
S 584 23 0 0 0 6 611 582 4680 4 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 ngl
S 585 23 0 0 0 6 612 582 4684 4 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 nlon
S 587 23 0 0 0 6 683 582 4699 4 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 nout
S 588 23 0 0 0 6 684 582 4704 4 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 nerr
S 589 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 18 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 590 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 591 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 592 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 17 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 593 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 594 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 611 6 15 mo_control ngl
R 612 6 16 mo_control nlon
R 683 16 2 mo_doctor nout
R 684 16 3 mo_doctor nerr
S 697 7 6 0 0 73 1 582 5332 10a00024 51 A 0 0 0 0 0 0 702 0 0 0 704 0 0 0 0 0 0 0 0 701 0 0 703 582 0 0 0 0 aland
S 698 6 4 0 0 6 699 582 5338 40800006 0 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 z_b_0_1
S 699 6 4 0 0 6 708 582 5346 40800006 0 A 0 0 0 0 0 4 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 z_b_1
S 700 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 0 0 0 0 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 701 8 4 0 0 76 717 582 5352 40822004 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 aland$sd
S 702 6 4 0 0 7 703 582 5361 40802001 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 aland$p
S 703 6 4 0 0 7 701 582 5369 40802000 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 aland$o
S 704 22 1 0 0 9 1 582 5377 40000000 1000 A 0 0 0 0 0 0 0 697 0 0 0 0 701 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 aland$arrdsc
S 705 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 0 0 0 0 0 19 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 706 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 0 0 0 0 0 20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 707 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 0 0 0 0 0 23 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 708 7 4 0 4 79 714 582 5390 800024 100 A 0 0 0 0 0 16 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 itopo
S 709 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 720 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 710 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 121 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 711 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 87120 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 712 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 0 0 0 0 0 721 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 713 7 6 0 0 82 1 582 5396 10a00024 51 A 0 0 0 0 0 0 717 0 0 0 719 0 0 0 0 0 0 0 0 716 0 0 718 582 0 0 0 0 iland
S 714 6 4 0 0 6 715 582 5402 40800006 0 A 0 0 0 0 0 348560 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 z_b_2
S 715 6 4 0 0 6 721 582 5408 40800006 0 A 0 0 0 0 0 348564 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 z_b_3
S 716 8 4 0 0 85 724 582 5414 40822004 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 iland$sd
S 717 6 4 0 0 7 718 582 5423 40802001 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 iland$p
S 718 6 4 0 0 7 716 582 5431 40802000 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 iland$o
S 719 22 1 0 0 6 1 582 5439 40000000 1000 A 0 0 0 0 0 0 0 713 0 0 0 0 716 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 iland$arrdsc
S 720 7 6 0 0 88 1 582 5452 10a00024 59 A 0 0 0 0 0 0 724 0 0 0 726 0 0 0 0 0 0 0 0 723 0 0 725 582 0 0 0 0 bzone
S 721 6 4 0 0 6 722 582 5458 40800006 0 A 0 0 0 0 0 348568 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 z_b_4
S 722 6 4 0 0 6 728 582 5464 40800006 0 A 0 0 0 0 0 348572 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 z_b_5
S 723 8 4 0 0 91 730 582 5470 40822004 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 bzone$sd
S 724 6 4 0 0 7 725 582 5479 40802001 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 bzone$p
S 725 6 4 0 0 7 723 582 5487 40802000 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 bzone$o
S 726 22 1 0 0 9 1 582 5495 40000000 1000 A 0 0 0 0 0 0 0 720 0 0 0 0 723 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 bzone$arrdsc
S 727 7 6 0 0 94 1 582 5508 10a00024 51 A 0 0 0 0 0 0 730 0 0 0 732 0 0 0 0 0 0 0 0 729 0 0 731 582 0 0 0 0 agx
S 728 6 4 0 0 6 734 582 5512 40800006 0 A 0 0 0 0 0 348576 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 z_b_6
S 729 8 4 0 0 97 736 582 5518 40822004 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 agx$sd
S 730 6 4 0 0 7 731 582 5525 40802001 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 agx$p
S 731 6 4 0 0 7 729 582 5531 40802000 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 agx$o
S 732 22 1 0 0 9 1 582 5537 40000000 1000 A 0 0 0 0 0 0 0 727 0 0 0 0 729 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 agx$arrdsc
S 733 7 6 0 0 100 1 582 5548 10a00024 51 A 0 0 0 0 0 0 736 0 0 0 738 0 0 0 0 0 0 0 0 735 0 0 737 582 0 0 0 0 agy
S 734 6 4 0 0 6 742 582 5552 40800006 0 A 0 0 0 0 0 348580 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 z_b_7
S 735 8 4 0 0 103 745 582 5558 40822004 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 agy$sd
S 736 6 4 0 0 7 737 582 5565 40802001 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 agy$p
S 737 6 4 0 0 7 735 582 5571 40802000 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 agy$o
S 738 22 1 0 0 9 1 582 5577 40000000 1000 A 0 0 0 0 0 0 0 733 0 0 0 0 735 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 agy$arrdsc
S 739 7 4 0 4 106 740 582 5588 800024 100 A 0 0 0 0 0 0 0 0 0 0 0 0 749 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 sgx
S 740 7 4 0 4 109 1 582 5592 800024 100 A 0 0 0 0 0 5760 0 0 0 0 0 0 749 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 sgy
S 741 7 6 0 0 112 1 582 5596 10a00024 59 A 0 0 0 0 0 0 745 0 0 0 747 0 0 0 0 0 0 0 0 744 0 0 746 582 0 0 0 0 sstatm
S 742 6 4 0 0 6 743 582 5603 40800006 0 A 0 0 0 0 0 348584 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 z_b_8
S 743 6 4 0 0 6 1 582 5609 40800006 0 A 0 0 0 0 0 348588 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 z_b_9
S 744 8 4 0 0 115 698 582 5615 40822004 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 sstatm$sd
S 745 6 4 0 0 7 746 582 5625 40802001 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 sstatm$p
S 746 6 4 0 0 7 744 582 5634 40802000 1020 A 0 0 0 0 0 0 0 0 0 0 0 0 748 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 sstatm$o
S 747 22 1 0 0 9 1 582 5643 40000000 1000 A 0 0 0 0 0 0 0 741 0 0 0 0 744 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 sstatm$arrdsc
S 748 11 0 0 4 9 696 582 5657 40800000 801000 A 0 0 0 0 0 349216 0 0 702 743 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _mo_landsea$0
S 749 11 0 0 4 9 748 582 5671 40800000 801000 A 0 0 0 0 0 6728 0 0 739 740 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _mo_landsea$2
S 750 23 5 0 0 0 754 582 5685 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 amask
S 751 7 3 0 0 124 1 750 5691 20000004 10003000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 slmm
S 752 7 3 0 0 118 1 750 5696 800204 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 vlat
S 753 7 3 0 0 121 1 750 5701 800204 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 vlon
S 754 14 5 0 0 0 1 750 5685 20000200 400000 A 0 0 0 0 0 0 0 2 3 0 0 0 0 0 0 0 0 0 0 0 0 20 0 582 0 0 0 0 amask
F 754 3 751 752 753
S 755 6 1 0 0 6 1 750 5706 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_124
S 756 6 1 0 0 6 1 750 5714 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_126
S 757 6 1 0 0 6 1 750 5722 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0_2
S 758 6 1 0 0 6 1 750 5730 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 759 6 1 0 0 6 1 750 5738 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3_1
S 760 6 1 0 0 6 1 750 5746 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5_1
S 761 6 1 0 0 6 1 750 5754 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6_1
S 762 6 1 0 0 6 1 750 5762 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_142
S 763 6 1 0 0 6 1 750 5770 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_145
S 764 23 5 0 0 0 773 582 5778 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 view
S 765 7 3 0 0 127 1 764 5783 800204 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 fland
S 766 1 3 0 0 6 1 764 5789 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 is
S 767 1 3 0 0 6 1 764 5792 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ie
S 768 1 3 0 0 6 1 764 5795 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 js
S 769 1 3 0 0 6 1 764 5798 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 je
S 770 1 3 0 0 6 1 764 5801 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 jv
S 771 6 3 0 0 6 1 764 4684 800004 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nlon
S 772 6 3 0 0 6 1 764 4680 800004 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ngl
S 773 14 5 0 0 0 1 764 5778 200 400000 A 0 0 0 0 0 0 0 6 8 0 0 0 0 0 0 0 0 0 0 0 0 432 0 582 0 0 0 0 view
F 773 8 765 771 772 766 767 768 769 770
S 774 6 1 0 0 6 1 764 5804 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_133
S 775 6 1 0 0 6 1 764 5812 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_134
S 776 6 1 0 0 6 1 764 5820 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_137
S 777 6 1 0 0 6 1 764 5828 40800006 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_139
S 778 23 5 0 0 0 785 582 5836 0 0 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 viewtop
S 779 7 3 0 0 130 1 778 5844 800004 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 top
S 780 1 3 0 0 6 1 778 5789 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 is
S 781 1 3 0 0 6 1 778 5792 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ie
S 782 1 3 0 0 6 1 778 5795 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 js
S 783 1 3 0 0 6 1 778 5798 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 je
S 784 1 3 0 0 6 1 778 5801 4 3000 A 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 jv
S 785 14 5 0 0 0 1 778 5836 0 400000 A 0 0 0 0 0 0 0 15 6 0 0 0 0 0 0 0 0 0 0 0 0 461 0 582 0 0 0 0 viewtop
F 785 6 779 780 781 782 783 784
A 12 2 0 0 0 6 589 0 0 0 12 0 0 0 0 0 0 0 0 0
A 14 2 0 0 0 6 594 0 0 0 14 0 0 0 0 0 0 0 0 0
A 16 2 0 0 0 6 590 0 0 0 16 0 0 0 0 0 0 0 0 0
A 18 2 0 0 0 6 591 0 0 0 18 0 0 0 0 0 0 0 0 0
A 22 2 0 0 0 6 592 0 0 0 22 0 0 0 0 0 0 0 0 0
A 24 2 0 0 0 6 593 0 0 0 24 0 0 0 0 0 0 0 0 0
A 38 2 0 0 0 6 700 0 0 0 38 0 0 0 0 0 0 0 0 0
A 39 1 0 3 0 76 701 0 0 0 0 0 0 0 0 0 0 0 0 0
A 40 10 0 0 0 6 39 4 0 0 0 0 0 0 0 0 0 0 0 0
X 1 16
A 41 10 0 0 40 6 39 7 0 0 0 0 0 0 0 0 0 0 0 0
X 1 18
A 42 4 0 0 0 6 41 0 3 0 0 0 0 2 0 0 0 0 0 0
A 43 4 0 0 0 6 40 0 42 0 0 0 0 1 0 0 0 0 0 0
A 44 2 0 0 0 6 705 0 0 0 44 0 0 0 0 0 0 0 0 0
A 45 10 0 0 41 6 39 16 0 0 0 0 0 0 0 0 0 0 0 0
X 1 44
A 46 2 0 0 0 6 706 0 0 0 46 0 0 0 0 0 0 0 0 0
A 47 10 0 0 45 6 39 19 0 0 0 0 0 0 0 0 0 0 0 0
X 1 46
A 48 4 0 0 0 6 47 0 3 0 0 0 0 2 0 0 0 0 0 0
A 49 4 0 0 0 6 45 0 48 0 0 0 0 1 0 0 0 0 0 0
A 50 10 0 0 47 6 39 10 0 0 0 0 0 0 0 0 0 0 0 0
X 1 22
A 51 2 0 0 0 6 707 0 0 0 51 0 0 0 0 0 0 0 0 0
A 52 10 0 0 50 6 39 22 0 0 0 0 0 0 0 0 0 0 0 0
X 1 51
A 53 10 0 0 52 6 39 13 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 54 10 0 0 53 6 39 1 0 0 0 0 0 0 0 0 0 0 0 0
X 1 14
A 55 2 0 0 0 6 709 0 0 0 55 0 0 0 0 0 0 0 0 0
A 56 2 0 0 0 6 710 0 0 0 56 0 0 0 0 0 0 0 0 0
A 57 2 0 0 0 6 711 0 0 0 57 0 0 0 0 0 0 0 0 0
A 58 2 0 0 0 6 712 0 0 0 58 0 0 0 0 0 0 0 0 0
A 61 1 0 3 0 85 716 0 0 0 0 0 0 0 0 0 0 0 0 0
A 62 10 0 0 0 6 61 4 0 0 0 0 0 0 0 0 0 0 0 0
X 1 16
A 63 10 0 0 62 6 61 7 0 0 0 0 0 0 0 0 0 0 0 0
X 1 18
A 64 4 0 0 0 6 63 0 3 0 0 0 0 2 0 0 0 0 0 0
A 65 4 0 0 0 6 62 0 64 0 0 0 0 1 0 0 0 0 0 0
A 66 10 0 0 63 6 61 16 0 0 0 0 0 0 0 0 0 0 0 0
X 1 44
A 67 10 0 0 66 6 61 19 0 0 0 0 0 0 0 0 0 0 0 0
X 1 46
A 68 4 0 0 0 6 67 0 3 0 0 0 0 2 0 0 0 0 0 0
A 69 4 0 0 27 6 66 0 68 0 0 0 0 1 0 0 0 0 0 0
A 70 10 0 0 67 6 61 10 0 0 0 0 0 0 0 0 0 0 0 0
X 1 22
A 71 10 0 0 70 6 61 22 0 0 0 0 0 0 0 0 0 0 0 0
X 1 51
A 72 10 0 0 71 6 61 13 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 73 10 0 0 72 6 61 1 0 0 0 0 0 0 0 0 0 0 0 0
X 1 14
A 76 1 0 3 0 91 723 0 0 0 0 0 0 0 0 0 0 0 0 0
A 77 10 0 0 0 6 76 4 0 0 0 0 0 0 0 0 0 0 0 0
X 1 16
A 78 10 0 0 77 6 76 7 0 0 0 0 0 0 0 0 0 0 0 0
X 1 18
A 79 4 0 0 0 6 78 0 3 0 0 0 0 2 0 0 0 0 0 0
A 80 4 0 0 0 6 77 0 79 0 0 0 0 1 0 0 0 0 0 0
A 81 10 0 0 78 6 76 16 0 0 0 0 0 0 0 0 0 0 0 0
X 1 44
A 82 10 0 0 81 6 76 19 0 0 0 0 0 0 0 0 0 0 0 0
X 1 46
A 83 4 0 0 30 6 82 0 3 0 0 0 0 2 0 0 0 0 0 0
A 84 4 0 0 0 6 81 0 83 0 0 0 0 1 0 0 0 0 0 0
A 85 10 0 0 82 6 76 10 0 0 0 0 0 0 0 0 0 0 0 0
X 1 22
A 86 10 0 0 85 6 76 22 0 0 0 0 0 0 0 0 0 0 0 0
X 1 51
A 87 10 0 0 86 6 76 13 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 88 10 0 0 87 6 76 1 0 0 0 0 0 0 0 0 0 0 0 0
X 1 14
A 90 1 0 1 0 97 729 0 0 0 0 0 0 0 0 0 0 0 0 0
A 91 10 0 0 0 6 90 4 0 0 0 0 0 0 0 0 0 0 0 0
X 1 16
A 92 10 0 0 91 6 90 7 0 0 0 0 0 0 0 0 0 0 0 0
X 1 18
A 93 4 0 0 0 6 92 0 3 0 0 0 0 2 0 0 0 0 0 0
A 94 4 0 0 0 6 91 0 93 0 0 0 0 1 0 0 0 0 0 0
A 95 10 0 0 92 6 90 10 0 0 0 0 0 0 0 0 0 0 0 0
X 1 22
A 96 10 0 0 95 6 90 13 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 97 10 0 0 96 6 90 1 0 0 0 0 0 0 0 0 0 0 0 0
X 1 14
A 99 1 0 1 0 103 735 0 0 0 0 0 0 0 0 0 0 0 0 0
A 100 10 0 0 0 6 99 4 0 0 0 0 0 0 0 0 0 0 0 0
X 1 16
A 101 10 0 0 100 6 99 7 0 0 0 0 0 0 0 0 0 0 0 0
X 1 18
A 102 4 0 0 16 6 101 0 3 0 0 0 0 2 0 0 0 0 0 0
A 103 4 0 0 0 6 100 0 102 0 0 0 0 1 0 0 0 0 0 0
A 104 10 0 0 101 6 99 10 0 0 0 0 0 0 0 0 0 0 0 0
X 1 22
A 105 10 0 0 104 6 99 13 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 106 10 0 0 105 6 99 1 0 0 0 0 0 0 0 0 0 0 0 0
X 1 14
A 109 1 0 3 0 115 744 0 0 0 0 0 0 0 0 0 0 0 0 0
A 110 10 0 0 0 6 109 4 0 0 0 0 0 0 0 0 0 0 0 0
X 1 16
A 111 10 0 0 110 6 109 7 0 0 0 0 0 0 0 0 0 0 0 0
X 1 18
A 112 4 0 0 0 6 111 0 3 0 0 0 0 2 0 0 0 0 0 0
A 113 4 0 0 0 6 110 0 112 0 0 0 0 1 0 0 0 0 0 0
A 114 10 0 0 111 6 109 16 0 0 0 0 0 0 0 0 0 0 0 0
X 1 44
A 115 10 0 0 114 6 109 19 0 0 0 0 0 0 0 0 0 0 0 0
X 1 46
A 116 4 0 0 0 6 115 0 3 0 0 0 0 2 0 0 0 0 0 0
A 117 4 0 0 0 6 114 0 116 0 0 0 0 1 0 0 0 0 0 0
A 118 10 0 0 115 6 109 10 0 0 0 0 0 0 0 0 0 0 0 0
X 1 22
A 119 10 0 0 118 6 109 22 0 0 0 0 0 0 0 0 0 0 0 0
X 1 51
A 120 10 0 0 119 6 109 13 0 0 0 0 0 0 0 0 0 0 0 0
X 1 24
A 121 10 0 0 120 6 109 1 0 0 0 0 0 0 0 0 0 0 0 0
X 1 14
A 122 1 0 0 0 6 611 0 0 0 0 0 0 0 0 0 0 0 0 0
A 123 1 0 0 0 6 755 0 0 0 0 0 0 0 0 0 0 0 0 0
A 124 1 0 0 0 6 612 0 0 0 0 0 0 0 0 0 0 0 0 0
A 125 1 0 0 0 6 756 0 0 0 0 0 0 0 0 0 0 0 0 0
A 126 1 0 0 0 6 761 0 0 0 0 0 0 0 0 0 0 0 0 0
A 127 1 0 0 0 6 757 0 0 0 0 0 0 0 0 0 0 0 0 0
A 128 1 0 0 0 6 762 0 0 0 0 0 0 0 0 0 0 0 0 0
A 129 1 0 0 0 6 759 0 0 0 0 0 0 0 0 0 0 0 0 0
A 130 1 0 0 0 6 758 0 0 0 0 0 0 0 0 0 0 0 0 0
A 131 1 0 0 0 6 763 0 0 0 0 0 0 0 0 0 0 0 0 0
A 132 1 0 0 0 6 760 0 0 0 0 0 0 0 0 0 0 0 0 0
A 133 1 0 0 0 6 777 0 0 0 0 0 0 0 0 0 0 0 0 0
A 134 1 0 0 0 6 771 0 0 0 0 0 0 0 0 0 0 0 0 0
A 135 1 0 0 0 6 774 0 0 0 0 0 0 0 0 0 0 0 0 0
A 136 1 0 0 0 6 772 0 0 0 0 0 0 0 0 0 0 0 0 0
A 137 1 0 0 0 6 775 0 0 0 0 0 0 0 0 0 0 0 0 0
A 138 1 0 0 0 6 776 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
Z