#include "gtest/gtest.h"

#include "GCP.h"
#include "utils_test.h"

using namespace GCP;

TEST(PercentileQueries, HugeTestWithNoDuplicates){
  test_several_confidences(
      {121958, 671155, 131932, 365838, 259178, 644167, 110268, 732180, 54886, 137337, 999890, 521430, 954698, 87498, 899159, 912756, 175203, 191335, 278167, 41090, 329365, 64820, 787201, 321879, 718315, 327069, 776997, 199041, 791743, 103355, 235796, 214176, 184779, 347449, 421909, 989436, 258795, 486232, 917040, 500186, 156730, 870910, 384681, 149503, 654811, 527035, 648143, 452366, 65725, 129981, 84654, 953277, 905778, 591723, 319030, 328947, 555839, 902648, 273538, 331236, 528178, 565894, 489492, 256840, 473254, 349457, 665987, 270936, 239931, 239629, 205041, 698361, 731912, 417113, 68148, 777089, 648531, 251995, 906606, 678843, 910790, 572843, 256508, 803591, 896942, 106530, 604365, 822352, 460337, 724839, 805889, 459773, 208261, 764469, 341097, 315139, 171829, 271836, 438974, 989913, 202283, 196769, 561353, 223165, 623587, 950110, 535822, 487879, 564685, 782038, 353531, 263160, 579879, 220884, 23247, 24300, 467281, 607086, 533556, 348951, 274329, 719064, 742139, 432315, 926075, 516588, 236584, 562332, 837646, 978732, 376896, 48984, 302918, 264712, 120151, 455808, 985067, 133767, 648663, 703550, 136330, 480754, 301648, 688519, 13986, 586146, 129312, 12666, 300804, 134633, 288998, 850937, 628776, 707611, 284806, 794824, 30535, 996107, 375713, 838688, 351279, 913910, 821654, 601661, 997079, 647972, 489570, 456551, 273109, 24538, 894498, 201664, 491234, 899684, 139182, 472525, 158338, 184064, 214020, 623094, 924414, 897421, 615270, 645914, 119176, 917710, 218126, 218969, 122409, 547707, 968140, 253618, 769598, 50015, 772838, 122096, 631347, 185340, 110687, 874371, 411357, 735716, 105878, 727270, 973548, 198286, 773290, 774684, 164899, 973324, 324767, 999238, 861882, 879858, 379989, 415515, 43585, 371369, 64044, 960061, 442296, 109556, 526981, 881691, 767595, 200235, 25939, 966429, 960445, 152906, 542335, 288249, 75766, 408923, 415192, 568550, 387261, 810208, 685440, 652664, 984346, 454589, 137848, 617075, 127948, 579304, 402690, 499046, 451269, 21959, 660890, 750216, 905533, 790180, 228576, 537833, 914091, 983703, 165838, 335674, 879989, 61087, 593128, 118451, 447600, 837437, 439792, 39353, 193075, 214283, 314877, 65318, 765313, 482690, 416880, 591460, 740976, 781239, 93264, 770746, 813168, 225281, 637147, 898613, 394070, 352228, 431839, 269536, 317824, 257426, 177789, 513153, 742452, 837291, 125657, 235167, 681669, 480671, 719094, 70467, 183734, 622794, 109751, 381050, 930192, 811774, 692517, 461079, 368452, 969587, 106081, 983237, 613333, 558986, 110078, 867055, 22671, 418400, 355528, 959611, 610490, 51663, 15708, 878338, 311955, 133883, 318394, 112547, 989586, 910681, 559042, 751102, 762002, 142483, 658271, 984774, 214835, 401896, 84896, 477095, 71295, 743974, 887121, 273255, 299648, 657162, 8155, 466872, 825816, 374705, 571542, 808350, 610269, 978217, 554594, 147718, 308987, 645263, 358745, 107512, 147443, 470587, 373616, 657409, 686976, 365871, 375037, 807364, 117796, 291999, 192506, 990198, 698376, 269544, 380002, 698002, 539439, 968911, 590978, 528787, 530583, 907317, 388215, 321184, 798615, 422515, 668234, 378480, 32711, 791971, 501157, 184423, 978771, 516349, 703714, 912495, 657917, 824792, 737378, 133272, 767836, 986001, 702335, 807789, 643409, 77505, 2869, 781474, 454351, 312252, 636584, 618467, 90272, 38467, 809760, 877844, 480047, 659347, 58871, 939903, 128391, 439430, 610704, 856703, 830496, 256687, 725451, 107450, 171890, 438741, 281974, 934933, 400109, 116381, 473125, 571621, 134508, 729650, 862645, 399111, 150810, 292890, 573665, 965908, 840477, 360032, 957294, 369599, 660960, 750020, 879611, 171836, 729903, 45714, 974339, 627234, 896048, 630271, 466960, 711851, 124123, 936093, 305628, 82989, 298356, 829957, 178274, 687995, 325352, 820260, 411927, 789852, 883184, 131373, 438452, 115294, 459451, 971744, 648307, 65726, 141564, 377812, 268246, 705696, 118015, 740674, 254079, 550929, 689944, 688105, 513758, 940597, 983609, 365890, 445101, 36631, 766577, 72991, 4014, 273237, 18070, 822209, 56958, 82074, 535017, 755073, 176089, 197392, 158823, 630152, 393002, 733429, 498863, 367654, 579036, 530089, 150262, 442905, 235362, 444209, 553880, 797079, 664076, 219963, 480761, 248710, 220984, 897059, 856748, 782355, 264512, 601863, 268799, 237549, 594319, 993933, 484299, 612210, 448982, 124046, 406619, 32097, 996161, 666326, 649150, 230229, 118834, 93848, 52921, 705086, 650429, 705660, 284821, 55609, 227897, 974165, 31024, 70313, 389957, 184078, 357429, 932842, 461243, 113632, 997639, 937012, 251451, 3051, 704107, 376836, 611523, 971525, 721772, 189407, 449395, 949597, 462894, 211810, 957238, 727975, 73523, 760463, 94476, 233841, 803451, 870045, 222866, 86672, 598135, 344894, 734994, 28251, 627769, 374710, 550233, 877284, 32217, 8308, 268093, 183062, 394366, 744840, 100235, 74740, 378496, 853049, 750201, 803328, 917230, 177247, 663165, 529525, 20056, 176615, 187628, 638836, 373632, 134415, 406716, 446438, 553663, 665022, 541252, 637717, 252764, 714998, 466882, 740427, 512153, 67215, 855474, 144356, 203861, 355612, 454605, 9435, 447556, 699438, 865264, 575581, 78781, 191475, 160196, 226507, 306063, 764121, 776665, 690607, 121172, 356902, 484963, 455968, 318717, 59101, 495972, 288790, 890617, 547337, 734965, 360260, 883811, 592673, 650419, 433374, 720777, 547577, 372370, 887225, 271967, 191232, 304119, 9540, 540652, 397827, 551951, 51991, 188751, 366326, 686337, 618992, 49115, 353151, 743583, 191962, 721747, 149143, 714635, 989873, 536226, 952699, 353824, 4000, 725948, 887474, 355754, 251492, 890123, 499522, 491072, 38304, 504999, 331593, 607274, 151836, 108940, 824926, 265517, 667521, 621820, 879090, 385777, 624674, 551638, 617296, 880729, 549639, 612956, 638873, 932809, 337497, 207869, 65953, 843890, 500328, 806790, 3267, 82745, 876020, 824945, 157164, 535626, 131484, 531831, 892323, 951000, 60692, 220156, 223992, 516771, 830857, 563044, 543176, 691095, 833727, 114352, 869474, 85999, 17955, 824273, 70640, 196582, 952415, 976535, 744086, 400573, 425695, 61476, 306955, 167280, 675510, 595468, 698646, 637144, 854882, 414568, 186141, 898320, 634736, 538685, 959059, 178031, 918232, 968021, 119180, 66234, 564242, 178352, 79459, 409995, 195004, 466152, 941970, 954443, 725256, 952774, 59163, 475853, 51934, 720255, 425523, 840530, 454137, 5486, 139407, 873332, 650308, 627083, 220184, 666490, 349652, 761315, 642612, 322710, 197775, 289336, 787494, 608364, 89780, 152617, 545977, 797606, 350605, 912606, 889465, 533636, 941218, 38102, 74460, 364778, 89930, 531089, 627438, 727627, 179208, 154697, 213945, 231915, 727952, 208124, 437477, 525830, 283821, 177804, 705703, 39081, 718216, 355249, 219930, 805497, 517313, 214532, 612380, 519332, 316837, 289106, 28295, 13807, 435564, 157504, 238067, 485008, 975358, 957766, 692440, 394540, 877699, 517923, 166981, 744699, 19870, 952082, 203196, 248683, 234677, 26790, 643162, 811721, 448345, 749074, 487462, 539773, 830914, 86188, 529036, 742129, 664795, 790141, 648249, 30355, 711131, 451015, 406332, 156542, 173416, 305062, 838656, 563586, 501708, 866651, 61629, 328184, 492734, 755195, 290160, 548248, 820023, 213920, 889765, 242804, 376965, 460857, 822827, 504492, 653983, 415916, 821052, 758574, 128148, 806991, 486261, 355540, 192714, 731899, 283876, 259747, 615650, 423570, 275987, 638353, 858158, 17640, 52528, 824845, 259214, 839070, 970240, 501108, 634421, 368501, 724226, 150159, 172502, 155576, 924618, 865291, 886089, 928479, 11023, 545735, 388043, 286103, 66203, 905588, 160775, 9337, 570715, 514723, 323673, 364167, 346809, 630075, 687281, 295451, 420571, 127016, 791267, 544447, 497562, 702965, 383934, 892048, 955592, 865184, 317428, 319187, 690907, 883439, 82844, 551820, 165421, 122402, 118012, 438661, 398929, 580466, 573695, 463556, 575278, 993176, 129473, 983049, 478007, 844573, 417009, 318714, 226156, 168964, 525303, 61813, 27712, 624785, 925919, 767214, 112816, 128778, 586580, 285977, 135230, 718040, 806741},
      {0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000, 23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000, 32000, 33000, 34000, 35000, 36000, 37000, 38000, 39000, 40000, 41000, 42000, 43000, 44000, 45000, 46000, 47000, 48000, 49000, 50000, 51000, 52000, 53000, 54000, 55000, 56000, 57000, 58000, 59000, 60000, 61000, 62000, 63000, 64000, 65000, 66000, 67000, 68000, 69000, 70000, 71000, 72000, 73000, 74000, 75000, 76000, 77000, 78000, 79000, 80000, 81000, 82000, 83000, 84000, 85000, 86000, 87000, 88000, 89000, 90000, 91000, 92000, 93000, 94000, 95000, 96000, 97000, 98000, 99000, 100000, 101000, 102000, 103000, 104000, 105000, 106000, 107000, 108000, 109000, 110000, 111000, 112000, 113000, 114000, 115000, 116000, 117000, 118000, 119000, 120000, 121000, 122000, 123000, 124000, 125000, 126000, 127000, 128000, 129000, 130000, 131000, 132000, 133000, 134000, 135000, 136000, 137000, 138000, 139000, 140000, 141000, 142000, 143000, 144000, 145000, 146000, 147000, 148000, 149000, 150000, 151000, 152000, 153000, 154000, 155000, 156000, 157000, 158000, 159000, 160000, 161000, 162000, 163000, 164000, 165000, 166000, 167000, 168000, 169000, 170000, 171000, 172000, 173000, 174000, 175000, 176000, 177000, 178000, 179000, 180000, 181000, 182000, 183000, 184000, 185000, 186000, 187000, 188000, 189000, 190000, 191000, 192000, 193000, 194000, 195000, 196000, 197000, 198000, 199000, 200000, 201000, 202000, 203000, 204000, 205000, 206000, 207000, 208000, 209000, 210000, 211000, 212000, 213000, 214000, 215000, 216000, 217000, 218000, 219000, 220000, 221000, 222000, 223000, 224000, 225000, 226000, 227000, 228000, 229000, 230000, 231000, 232000, 233000, 234000, 235000, 236000, 237000, 238000, 239000, 240000, 241000, 242000, 243000, 244000, 245000, 246000, 247000, 248000, 249000, 250000, 251000, 252000, 253000, 254000, 255000, 256000, 257000, 258000, 259000, 260000, 261000, 262000, 263000, 264000, 265000, 266000, 267000, 268000, 269000, 270000, 271000, 272000, 273000, 274000, 275000, 276000, 277000, 278000, 279000, 280000, 281000, 282000, 283000, 284000, 285000, 286000, 287000, 288000, 289000, 290000, 291000, 292000, 293000, 294000, 295000, 296000, 297000, 298000, 299000, 300000, 301000, 302000, 303000, 304000, 305000, 306000, 307000, 308000, 309000, 310000, 311000, 312000, 313000, 314000, 315000, 316000, 317000, 318000, 319000, 320000, 321000, 322000, 323000, 324000, 325000, 326000, 327000, 328000, 329000, 330000, 331000, 332000, 333000, 334000, 335000, 336000, 337000, 338000, 339000, 340000, 341000, 342000, 343000, 344000, 345000, 346000, 347000, 348000, 349000, 350000, 351000, 352000, 353000, 354000, 355000, 356000, 357000, 358000, 359000, 360000, 361000, 362000, 363000, 364000, 365000, 366000, 367000, 368000, 369000, 370000, 371000, 372000, 373000, 374000, 375000, 376000, 377000, 378000, 379000, 380000, 381000, 382000, 383000, 384000, 385000, 386000, 387000, 388000, 389000, 390000, 391000, 392000, 393000, 394000, 395000, 396000, 397000, 398000, 399000, 400000, 401000, 402000, 403000, 404000, 405000, 406000, 407000, 408000, 409000, 410000, 411000, 412000, 413000, 414000, 415000, 416000, 417000, 418000, 419000, 420000, 421000, 422000, 423000, 424000, 425000, 426000, 427000, 428000, 429000, 430000, 431000, 432000, 433000, 434000, 435000, 436000, 437000, 438000, 439000, 440000, 441000, 442000, 443000, 444000, 445000, 446000, 447000, 448000, 449000, 450000, 451000, 452000, 453000, 454000, 455000, 456000, 457000, 458000, 459000, 460000, 461000, 462000, 463000, 464000, 465000, 466000, 467000, 468000, 469000, 470000, 471000, 472000, 473000, 474000, 475000, 476000, 477000, 478000, 479000, 480000, 481000, 482000, 483000, 484000, 485000, 486000, 487000, 488000, 489000, 490000, 491000, 492000, 493000, 494000, 495000, 496000, 497000, 498000, 499000, 500000, 501000, 502000, 503000, 504000, 505000, 506000, 507000, 508000, 509000, 510000, 511000, 512000, 513000, 514000, 515000, 516000, 517000, 518000, 519000, 520000, 521000, 522000, 523000, 524000, 525000, 526000, 527000, 528000, 529000, 530000, 531000, 532000, 533000, 534000, 535000, 536000, 537000, 538000, 539000, 540000, 541000, 542000, 543000, 544000, 545000, 546000, 547000, 548000, 549000, 550000, 551000, 552000, 553000, 554000, 555000, 556000, 557000, 558000, 559000, 560000, 561000, 562000, 563000, 564000, 565000, 566000, 567000, 568000, 569000, 570000, 571000, 572000, 573000, 574000, 575000, 576000, 577000, 578000, 579000, 580000, 581000, 582000, 583000, 584000, 585000, 586000, 587000, 588000, 589000, 590000, 591000, 592000, 593000, 594000, 595000, 596000, 597000, 598000, 599000, 600000, 601000, 602000, 603000, 604000, 605000, 606000, 607000, 608000, 609000, 610000, 611000, 612000, 613000, 614000, 615000, 616000, 617000, 618000, 619000, 620000, 621000, 622000, 623000, 624000, 625000, 626000, 627000, 628000, 629000, 630000, 631000, 632000, 633000, 634000, 635000, 636000, 637000, 638000, 639000, 640000, 641000, 642000, 643000, 644000, 645000, 646000, 647000, 648000, 649000, 650000, 651000, 652000, 653000, 654000, 655000, 656000, 657000, 658000, 659000, 660000, 661000, 662000, 663000, 664000, 665000, 666000, 667000, 668000, 669000, 670000, 671000, 672000, 673000, 674000, 675000, 676000, 677000, 678000, 679000, 680000, 681000, 682000, 683000, 684000, 685000, 686000, 687000, 688000, 689000, 690000, 691000, 692000, 693000, 694000, 695000, 696000, 697000, 698000, 699000, 700000, 701000, 702000, 703000, 704000, 705000, 706000, 707000, 708000, 709000, 710000, 711000, 712000, 713000, 714000, 715000, 716000, 717000, 718000, 719000, 720000, 721000, 722000, 723000, 724000, 725000, 726000, 727000, 728000, 729000, 730000, 731000, 732000, 733000, 734000, 735000, 736000, 737000, 738000, 739000, 740000, 741000, 742000, 743000, 744000, 745000, 746000, 747000, 748000, 749000, 750000, 751000, 752000, 753000, 754000, 755000, 756000, 757000, 758000, 759000, 760000, 761000, 762000, 763000, 764000, 765000, 766000, 767000, 768000, 769000, 770000, 771000, 772000, 773000, 774000, 775000, 776000, 777000, 778000, 779000, 780000, 781000, 782000, 783000, 784000, 785000, 786000, 787000, 788000, 789000, 790000, 791000, 792000, 793000, 794000, 795000, 796000, 797000, 798000, 799000, 800000, 801000, 802000, 803000, 804000, 805000, 806000, 807000, 808000, 809000, 810000, 811000, 812000, 813000, 814000, 815000, 816000, 817000, 818000, 819000, 820000, 821000, 822000, 823000, 824000, 825000, 826000, 827000, 828000, 829000, 830000, 831000, 832000, 833000, 834000, 835000, 836000, 837000, 838000, 839000, 840000, 841000, 842000, 843000, 844000, 845000, 846000, 847000, 848000, 849000, 850000, 851000, 852000, 853000, 854000, 855000, 856000, 857000, 858000, 859000, 860000, 861000, 862000, 863000, 864000, 865000, 866000, 867000, 868000, 869000, 870000, 871000, 872000, 873000, 874000, 875000, 876000, 877000, 878000, 879000, 880000, 881000, 882000, 883000, 884000, 885000, 886000, 887000, 888000, 889000, 890000, 891000, 892000, 893000, 894000, 895000, 896000, 897000, 898000, 899000, 900000, 901000, 902000, 903000, 904000, 905000, 906000, 907000, 908000, 909000, 910000, 911000, 912000, 913000, 914000, 915000, 916000, 917000, 918000, 919000, 920000, 921000, 922000, 923000, 924000, 925000, 926000, 927000, 928000, 929000, 930000, 931000, 932000, 933000, 934000, 935000, 936000, 937000, 938000, 939000, 940000, 941000, 942000, 943000, 944000, 945000, 946000, 947000, 948000, 949000, 950000, 951000, 952000, 953000, 954000, 955000, 956000, 957000, 958000, 959000, 960000, 961000, 962000, 963000, 964000, 965000, 966000, 967000, 968000, 969000, 970000, 971000, 972000, 973000, 974000, 975000, 976000, 977000, 978000, 979000, 980000, 981000, 982000, 983000, 984000, 985000, 986000, 987000, 988000, 989000, 990000, 991000, 992000, 993000, 994000, 995000, 996000, 997000, 998000, 999000, 1000000},
      {0.0, 0.0, 0.0, 0.17197802197802198, 0.4, 0.5669836956521739, 0.6192581491195204, 0.6567253653053577, 0.6941925814911951, 0.8672497570456754, 1.131018206338503, 1.1984490896830748, 1.2594643944004869, 1.3292725679228747, 1.5008130081300812, 1.558885017421603, 1.6151138716356108, 1.6668737060041408, 1.8391304347826087, 1.9516666666666667, 2.0698924731182795, 2.149605885444036, 2.2057584269662924, 2.3571180555555555, 2.4715099715099713, 2.632976445396146, 2.70716803760282, 2.8227765726681127, 2.9534322820037104, 3.134223300970874, 3.1827669902912623, 3.3950920245398772, 3.490959925442684, 3.707372448979592, 3.7328826530612247, 3.7583928571428573, 3.78390306122449, 3.825084976206662, 3.8930659415363698, 4.186807817589576, 4.337248128957974, 4.394818652849741, 4.436472945891784, 4.476553106212425, 4.51949271958666, 4.566463128229215, 4.608746177370031, 4.639327217125382, 4.669908256880734, 4.712213740458015, 4.898333333333333, 4.959769417475728, 5.201675977653632, 5.404020356234097, 5.454910941475827, 5.515767634854772, 5.628984432913269, 5.702195504443283, 5.75446941975954, 5.8560869565217395, 6.0547416612164815, 6.177974683544304, 6.508381891528463, 6.5532048408785295, 6.598027790228596, 6.736144578313253, 7.118799999999999, 7.378083588175332, 7.484137191854233, 7.539353348729792, 7.58554272517321, 7.854961832061069, 7.941568396226415, 8.001691729323309, 8.150907150480256, 8.32534113060429, 8.41345600920069, 8.470960322024151, 8.538793103448276, 8.632300884955752, 8.720688336520077, 8.758929254302103, 8.797170172084131, 9.10066066066066, 9.16072072072072, 9.309428830462377, 9.4005291005291, 9.639709443099273, 9.721998247151621, 9.765819456617002, 9.92046783625731, 10.024331550802138, 10.057754010695186, 10.091176470588234, 10.224203821656051, 10.309098801875326, 10.326462927591596, 10.343827053307866, 10.361191179024136, 10.378555304740406, 10.395919430456678, 10.424519230769231, 10.456570512820512, 10.488621794871795, 10.525564803804993, 10.565200158541419, 10.660098522167488, 10.85108695652174, 11.034173669467787, 11.10974025974026, 11.376146788990827, 11.616827956989248, 11.670591397849462, 11.822549019607845, 11.95111111111111, 12.068789808917197, 12.164949402023918, 12.24374558303887, 12.394444444444446, 12.748538011695906, 12.984449021627189, 13.083153770812928, 13.230434782608695, 13.53610262675626, 13.597189981673793, 13.757170795306388, 13.8252391464312, 13.898822663723326, 14.026, 14.341573033707865, 14.601364942528734, 14.673204022988505, 14.905074626865671, 14.979701492537313, 15.221992481203007, 15.561474036850921, 15.67, 15.76653426017875, 15.911394302848576, 15.986356821589206, 16.127491886879927, 16.17385257301808, 16.247442872687703, 16.327602776294714, 16.38099305926321, 16.420861678004535, 16.453255587949464, 16.485649497894396, 16.61978947368421, 16.6899649122807, 16.87576219512195, 17.11851851851852, 17.22099871959027, 17.405248464544947, 17.46108319374651, 17.534470989761093, 17.643892339544514, 17.86221198156682, 18.059472422062353, 18.2128914785142, 18.28572469045885, 18.405455868089234, 18.42970417070805, 18.453952473326865, 18.478200775945684, 18.519348659003832, 18.714173228346457, 18.806354515050167, 18.942755344418053, 19.001256544502617, 19.036160558464225, 19.07106457242583, 19.31797385620915, 19.45448577680525, 19.532680470061557, 19.58864017907107, 19.689954853273136, 19.860917721518987, 20.186343612334802, 20.47570093457944, 20.520550077841204, 20.54649714582252, 20.57244421380384, 20.598391281785158, 20.78060606060606, 21.13939393939394, 21.28239700374532, 21.357767316745125, 21.43312555654497, 21.53795731707317, 21.632493150684933, 21.687287671232877, 22.006985294117648, 22.279224376731303, 22.347952306894765, 22.399792638672885, 22.46311787072243, 22.637078651685393, 22.8440313111546, 22.994569536423842, 23.08031825795645, 23.153533939818054, 23.25428109854604, 23.37853231106243, 23.51177966101695, 23.596525423728814, 23.63391089108911, 23.669271570014143, 23.751372549019607, 23.920822766976613, 23.94899971823049, 23.97717666948436, 24.009004739336493, 24.056398104265405, 24.273333333333333, 24.705013673655422, 24.735399574597388, 24.76578547553935, 24.796171376481315, 24.903225806451612, 25.11917098445596, 25.500850159404887, 25.553985122210413, 25.644816053511708, 25.800620636152058, 25.87820015515904, 25.98217142857143, 26.13546762589928, 26.21516936671576, 26.325650332728372, 26.386146400483966, 26.445729537366546, 26.50441329179647, 26.556334371754932, 26.61901913875598, 26.765918367346938, 27.025888324873097, 27.143108808290155, 27.287065637065638, 27.359731113956464, 27.502401670727462, 27.5372084928646, 27.572015315001742, 27.603333900323186, 27.620343595849636, 27.637353291376087, 27.654362986902537, 27.671372682428988, 27.688382377955435, 27.810580080262678, 27.84706311565122, 27.883546151039766, 28.10065019505852, 28.227634660421547, 28.38286334056399, 28.437916838205023, 28.47908604363936, 28.72730375426621, 28.84192841490139, 28.953524804177544, 29.2074128332845, 29.236712569586874, 29.266012305889248, 29.29531204219162, 29.362130177514793, 29.53577639751553, 29.618750000000002, 29.6575698757764, 29.696389751552793, 29.927272727272726, 30.132758620689657, 30.20711111111111, 30.402889667250438, 30.49045534150613, 30.858407079646017, 30.940470446320866, 31.00059633027523, 31.046467889908257, 31.092339449541285, 31.121880745994222, 31.14814814814815, 31.174415550302076, 31.201407688142933, 31.25554953979426, 31.41333333333333, 31.615484429065745, 31.71825396825397, 31.841798695246972, 31.888397017707362, 32.101851851851855, 32.380582524271844, 32.445676998368675, 32.50011223344557, 32.6042951971886, 32.64334244435767, 32.68238969152675, 32.71889845094665, 32.75332185886403, 32.78774526678141, 32.84984520123839, 32.930449826989616, 33.023222748815165, 33.12771653543307, 33.20682764363031, 33.29009159034138, 33.39342523860021, 33.585517241379314, 33.702214566929136, 33.75142716535433, 33.80043800539083, 33.8341307277628, 33.867823450134765, 33.915151515151514, 34.02849523809524, 34.06659047619048, 34.146946564885496, 34.25070671378092, 34.327580372250424, 34.53087719298246, 34.89041533546326, 35.0407110665999, 35.0907861792689, 35.21456077015644, 35.330114226375905, 35.42989031078611, 35.53982905982906, 35.6377402446127, 35.69598136284217, 35.783497757847535, 35.91267942583732, 36.03393907001603, 36.08738642437199, 36.2099730458221, 36.234476843910805, 36.25898064199951, 36.28348444008821, 36.317882611080634, 36.37273724629731, 36.41397222222222, 36.44175, 36.46952777777778, 36.497305555555556, 36.523781933105084, 36.550118514616805, 36.57645509612853, 36.605535248041775, 36.6577545691906, 36.72984375, 36.8366844207723, 36.9096837944664, 37.136516264428124, 37.25860534124629, 37.375974710221286, 37.48364030335861, 37.71235087719298, 37.78252631578947, 38.221428571428575, 38.31859582542694, 38.44338905775076, 38.51981351981352, 38.5975135975136, 38.71894036345022, 38.744535449193755, 38.77013053493729, 38.79572562068083, 38.92094339622641, 39.22522935779816, 39.35075301204819, 39.443358395989975, 39.64544626593807, 39.72265536723164, 39.77915254237288, 39.86303696303696, 39.95056179775281, 40.13429636533085, 40.38868501529052, 40.525556544968836, 40.80413223140496, 40.92814371257485, 41.1337575351641, 41.284615384615385, 41.39522900763359, 41.43294036061026, 41.46761442441054, 41.50883534136546, 41.62910583941606, 41.7150269541779, 41.78241239892183, 41.894501278772374, 42.045063145809415, 42.10141215106732, 42.13425287356322, 42.16709359605912, 42.19993431855501, 42.29344569288389, 42.51399452388196, 42.54441740188622, 42.57484027989048, 42.615698729582576, 42.73901098901099, 42.88907815631262, 43.03227513227513, 43.11309823677582, 43.20851180669962, 43.23596924766612, 43.26342668863262, 43.29088412959912, 43.512868146805616, 43.55817852288174, 43.6071828358209, 43.700367107195305, 43.7737885462555, 43.90276410450587, 43.94062854979175, 43.97849299507762, 44.06923076923077, 44.30871369294606, 44.49302325581395, 44.66891996891997, 44.72763703362506, 44.77369875633349, 44.832062780269055, 44.915016501650165, 45.04597156398104, 45.122017409114186, 45.17322068612391, 45.30496419270833, 45.321240234375, 45.33751627604166, 45.353792317708326, 45.370068359375, 45.38634440104166, 45.43382352941176, 45.564683663833804, 45.62858447488585, 45.674246575342465, 45.72279142707789, 45.77506534239414, 45.853641025641025, 46.20570175438597, 46.408306709265176, 46.4482428115016, 46.48817891373802, 46.60728527607362, 46.68397239263804, 46.78867713004484, 46.867240089753174, 46.95026833631485, 47.15369127516779, 47.30435835351089, 47.437345679012346, 47.49907407407407, 47.66663628076573, 47.735798983625074, 47.79226425748165, 48.13283458021613, 48.305488850771866, 48.41548275862069, 48.44996551724138, 48.48444827586207, 48.640248226950355, 48.86441441441441, 49.045850999394304, 49.11601208459215, 49.2171032357473, 49.25562403697997, 49.29414483821263, 49.61246105919003, 49.72174833635814, 49.751996370235936, 49.78224440411373, 49.82131062951496, 49.87291021671827, 49.979166666666664, 50.12870334744132, 50.16717968449404, 50.211835748792275, 50.29235104669887, 50.399232456140346, 50.44867647058823, 50.497696078431375, 50.812389839294966, 50.86423017107309, 50.919266625233064, 50.98141702921069, 51.18222222222222, 51.28104575163398, 51.46153205661948, 51.607501549907006, 51.6694978301302, 51.828628495339544, 51.89520639147803, 52.05106666666667, 52.10821494749846, 52.139098208770854, 52.16998147004324, 52.20176100628931, 52.264654088050314, 52.333666410453496, 52.47486338797814, 52.67198795180723, 52.886153846153846, 53.11048850574713, 53.14640804597701, 53.182327586206895, 53.300013978194016, 53.31399217221135, 53.32797036622868, 53.34194856024602, 53.35592675426335, 53.369904948280684, 53.38388314229802, 53.397861336315344, 53.4847, 53.62507772020726, 53.71703567035671, 53.778536285362854, 54.042250922509226, 54.205464868701206, 54.27643718949609, 54.33183984747378, 54.379504289799804, 54.41471727343145, 54.44053705138136, 54.46635682933127, 54.49217660728118, 54.61476976542137, 54.73518518518519, 54.88442694663167, 55.0855421686747, 55.28421985815603, 55.482411067193674, 55.609797101449274, 55.66776811594203, 55.82635771180304, 55.89876900796524, 56.144059405940595, 56.248164281269446, 56.319600938967135, 56.44177718832891, 56.625824800910124, 56.758, 56.869067405355494, 56.97907253269917, 57.0648308418568, 57.1429347826087, 57.301691176470584, 57.375220588235294, 57.654158964879855, 57.75406182602444, 57.86077441077441, 58.01001410437235, 58.30286214953271, 58.36127336448598, 58.516806722689076, 58.63261044176707, 58.70511598347633, 58.736892278360344, 58.76866857324436, 58.824999999999996, 58.94145391605365, 58.98472522717438, 59.06608784473953, 59.193820224719104, 59.36310975609756, 59.52605459057072, 59.60399096385542, 59.641641566265065, 59.6792921686747, 59.72078521939954, 59.76697459584295, 59.83446191051995, 60.031014729950904, 60.119099756691, 60.319267214150344, 60.38243840808591, 60.512127351664255, 60.541070911722144, 60.570014471780034, 60.59895803183792, 60.8206132879046, 60.909401408450705, 60.92700704225352, 60.94461267605634, 60.96221830985915, 60.97982394366197, 60.997429577464786, 61.10954979536153, 61.13228740336517, 61.15502501136881, 61.177762619372444, 61.2045643153527, 61.42915789473684, 61.571868131868136, 61.67321578505458, 61.75926892950392, 61.819947506561675, 61.85744281964754, 61.89493813273341, 61.92453204764605, 61.952892796369824, 61.98125354509359, 62.105475619504396, 62.14544364508394, 62.18541167066347, 62.22333700845277, 62.26008820286659, 62.29683939728041, 62.46660550458716, 62.53338582677166, 62.58587926509187, 62.83614163614163, 62.96943231441048, 63.21167108753316, 63.334434692823955, 63.386060918946825, 63.52456140350877, 63.59473684210526, 63.7601195559351, 63.9002828854314, 63.93564356435644, 63.97100424328147, 64.01848049281314, 64.16866666666667, 64.33799448022079, 64.50935596170584, 64.55287206266318, 64.59638816362053, 64.92293942403178, 65.0172440338722, 65.09422632794457, 65.36775092936803, 65.42124268054653, 65.45377358490566, 65.4863044892648, 65.61428571428571, 65.66839826839826, 65.77428571428571, 66.04449685534591, 66.303396630115, 66.33014174913079, 66.35688686814657, 66.38363198716235, 66.47054545454546, 66.67796833773087, 66.77600364963503, 66.9041788143829, 66.95276967930029, 67.01637426900585, 67.56919917864477, 67.67340241796201, 67.925548098434, 67.9702908277405, 68.02547384382108, 68.10205314009661, 68.20803913228414, 68.25057422373459, 68.29310931518502, 68.52344632768362, 68.66775092936803, 68.74232015554115, 68.90181405895692, 68.94716553287982, 68.99251700680273, 69.09165751920966, 69.29030837004404, 69.40383480825959, 69.64946653734239, 69.7671809256662, 69.82622389592605, 69.86045874700446, 69.89469359808285, 69.91940298507463, 69.94236509758898, 69.96532721010333, 69.98828932261767, 70.01470147014702, 70.04470447044704, 70.07470747074707, 70.10555555555555, 70.14094125973106, 70.17632696390658, 70.20877751259613, 70.23529567753911, 70.2618138424821, 70.28833200742508, 70.36243032329989, 70.50786885245901, 70.70454545454545, 70.93375438596492, 71.00844645550528, 71.24946808510639, 71.36728624535316, 71.50880583409298, 71.52703737465816, 71.54526891522335, 71.56350045578851, 71.5817319963537, 71.59996353691886, 71.94469696969698, 72.01939937866759, 72.05391784604764, 72.08843631342768, 72.2059829059829, 72.47277353689567, 72.59121552604698, 72.91556603773586, 72.96797693920335, 73.01105113636363, 73.03946022727273, 73.06786931818182, 73.0962784090909, 73.20535201149426, 73.24127155172414, 73.27719109195402, 73.4000657462196, 73.43293885601578, 73.46581196581197, 73.49868507560815, 73.79145527369826, 73.9780361757106, 74.12298969072165, 74.30929095354523, 74.35004074979625, 74.39079054604727, 74.53860911270984, 74.80393343419061, 74.87957639939486, 75.20149253731344, 75.26119402985074, 75.40485971943889, 75.45495991983968, 75.63283582089552, 75.76565252201762, 75.83717447916666, 76.00083102493075, 76.11708784596871, 76.17725631768953, 76.22040013119056, 76.25319776976058, 76.28599540833059, 76.50208152645273, 76.58881179531656, 76.84845269672856, 77.02321428571429, 77.30377893245158, 77.32739726027397, 77.35101558809636, 77.37463391591875, 77.39825224374115, 77.49788583509513, 77.78848758465011, 77.82261395114581, 77.84779652480483, 77.87297909846387, 77.8981616721229, 78.02382361645456, 78.05341817105652, 78.08301272565848, 78.1225516146109, 78.17548967707782, 78.26302816901408, 78.39970887918487, 78.44709768758848, 78.49428975932044, 78.66291469194313, 78.75435126582279, 78.86640502354788, 79.10930760499431, 79.1660612939841, 79.23501742160279, 79.31214149139579, 79.35994263862332, 79.4358407079646, 79.55093256814921, 79.61595153962645, 79.66643109540637, 79.80326086956522, 79.92195180722892, 79.94604819277109, 79.97014457831325, 79.99424096385542, 80.19326241134752, 80.3133099463475, 80.33394552208007, 80.35458109781264, 80.3752166735452, 80.39585224927777, 80.5214588634436, 80.56386768447837, 80.65121107266435, 80.8754369825207, 81.10101647388713, 81.13606729758149, 81.17111812127585, 81.20780487804878, 81.25215077605321, 81.29649667405765, 81.4390485629336, 81.508168894547, 81.52938680246127, 81.55060471037555, 81.57182261828983, 81.5930405262041, 81.82145855194123, 81.87392444910809, 82.01302816901408, 82.3024128686327, 82.53761140819964, 82.64609929078014, 82.75357142857143, 82.85234633179114, 83.01621233859397, 83.08794835007173, 83.11213712618526, 83.12672501823486, 83.14131291028447, 83.15590080233407, 83.17048869438366, 83.18507658643327, 83.19966447848286, 83.39343434343435, 83.56234234234233, 83.81196403872752, 83.88112033195021, 84.30631458094145, 84.4044433711664, 84.42859212750544, 84.45274088384448, 84.47688964018353, 84.50797773654917, 84.80305723426946, 84.8386064699609, 84.87415570565233, 84.90765993265994, 84.93571829405163, 84.96377665544333, 84.99183501683501, 85.2350495049505, 85.48167539267016, 85.56609808102345, 85.7139880952381, 85.74374999999999, 85.7735119047619, 85.81610541727672, 85.90670961659335, 85.92242300439976, 85.93813639220618, 85.95384978001258, 85.96956316781899, 85.98527655562539, 86.00298295454546, 86.05033143939394, 86.09767992424241, 86.1518821603928, 86.21993243243243, 86.34279902359643, 86.51787234042553, 86.58879432624113, 86.62261009667024, 86.64946294307197, 86.67631578947369, 86.71546526867628, 86.81398188263096, 86.8533674675069, 86.89275305238283, 87.15213235294118, 87.28638613861386, 87.33906572964035, 87.38040512608517, 87.49211908931699, 87.60371593724194, 87.64500412881915, 87.68629232039636, 87.76429258902792, 87.8381443298969, 87.89878714372347, 87.97753164556963, 88.23157894736842, 88.38803191489362, 88.70148648648649, 88.82817047817048, 88.92069658405894, 88.98767582049565, 89.20829675153644, 89.25219490781387, 89.29609306409131, 89.38827519379845, 89.62641888498241, 89.67664490205927, 89.86564245810057, 90.0267645003494, 90.09664570230608, 90.23112643678161, 90.27710344827587, 90.3323870967742, 90.39690322580645, 90.54957264957265, 90.76440489432703, 90.97087912087912, 91.110661268556, 91.14439946018894, 91.17813765182187, 91.21220103986136, 91.24686308492201, 91.28152512998267, 91.5268115942029, 91.65541490857946, 91.72030321046374, 91.75002972651606, 91.77975624256837, 91.91231671554253, 91.97096774193548, 92.221143847487, 92.34972375690607, 92.43082400813836, 92.46473380807053, 92.49864360800271, 92.75555555555556, 92.81242316402458, 92.82859915884826, 92.84477515367195, 92.86095114849563, 92.87712714331931, 92.893303138143, 93.0293620292083, 93.15192307692307, 93.2384775374376, 93.28007487520799, 93.33041447752481, 93.38879159369527, 93.43087504776462, 93.46908674054261, 93.60755619320899, 93.65538020086083, 93.70577586206896, 93.79198275862069, 93.89869423286181, 93.93417502594258, 93.96876513317191, 94.01397694524496, 94.16489533011273, 94.30039333945194, 94.31350465451685, 94.32661596958175, 94.33972728464666, 94.35283859971156, 94.36594991477646, 94.37906122984137, 94.39217254490626, 94.4785575048733, 94.6, 94.6924214417745, 95.04493041749502, 95.16200686106346, 95.3337807606264, 95.42478736330499, 95.48554070473877, 95.71809744779583, 95.7954369682908, 95.98644444444444, 96.1101592531576, 96.12846421380193, 96.14676917444628, 96.1650741350906, 96.18337909573495, 96.21765834932822, 96.33586683417086, 96.39868090452262, 96.61316568047337, 96.76324655436447, 96.85914396887159, 97.01620253164558, 97.07949367088607, 97.27325769854133, 97.4648675171737, 97.55454545454545, 97.62764565992865, 97.68709869203329, 97.90535296867695, 97.92872837774661, 97.95210378681627, 97.97547919588592, 97.99885460495558, 98.3461897356143, 98.57713310580205, 98.69989293361884, 98.72908296943231, 98.75819505094614, 98.78730713245997, 99.13052631578947, 99.22693082605775, 99.26051040967091, 99.29408999328408, 99.40308187672494, 99.44908003679853, 99.49507819687213, 99.69139433551199, 99.82257661038149, 99.88511569731082, 100.0}
  );
}


TEST(PercentileQueries, HugeTestWithDuplicates){
  test_several_confidences(
      {21, 1, 1, 11, 21, 26, 41, 46, 36, 6, 1, 46, 6, 46, 16, 1, 41, 46, 6, 16, 41, 6, 11, 36, 21, 26, 36, 46, 26, 41, 31, 31, 36, 31, 6, 21, 1, 41, 26, 26, 36, 41, 31, 31, 11, 1, 36, 16, 1, 1, 36, 31, 1, 31, 31, 16, 21, 31, 6, 41, 6, 1, 11, 36, 31, 1, 6, 16, 11, 21, 21, 26, 26, 46, 16, 36, 31, 41, 31, 46, 21, 1, 26, 41, 1, 16, 41, 6, 11, 1, 31, 21, 21, 21, 41, 36, 31, 46, 41, 16, 36, 11, 46, 31, 21, 11, 36, 21, 46, 11, 36, 31, 1, 31, 26, 11, 36, 11, 11, 1, 16, 6, 36, 31, 31, 36, 11, 16, 21, 11, 11, 6, 46, 11, 41, 46, 11, 16, 6, 41, 26, 1, 36, 21, 16, 6, 26, 31, 46, 36, 46, 11, 31, 11, 36, 26, 16, 46, 1, 31, 46, 21, 26, 16, 36, 36, 16, 31, 11, 16, 46, 36, 1, 26, 41, 16, 16, 16, 1, 46, 31, 26, 31, 41, 21, 11, 26, 46, 21, 46, 16, 46, 21, 1, 6, 46, 46, 36, 16, 21, 16, 46, 46, 21, 46, 41, 36, 36},
      {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48},
      {0.0, 5.288461538461538, 6.971153846153846, 8.653846153846153, 10.336538461538462, 12.01923076923077, 13.701923076923077, 15.384615384615385, 17.067307692307693, 18.75, 20.432692307692307, 22.115384615384617, 24.18269230769231, 26.250000000000004, 28.317307692307697, 30.384615384615387, 32.45192307692308, 34.519230769230774, 36.58653846153847, 38.65384615384616, 40.721153846153854, 42.78846153846154, 44.56730769230769, 46.34615384615385, 48.125, 49.90384615384615, 51.68269230769231, 53.65384615384615, 55.625, 57.59615384615385, 59.56730769230769, 61.53846153846154, 63.94230769230769, 66.34615384615384, 68.75, 71.15384615384616, 73.5576923076923, 75.57692307692308, 77.59615384615385, 79.61538461538463, 81.63461538461539, 83.65384615384616, 85.72115384615385, 87.78846153846155, 89.85576923076924, 91.92307692307693, 93.99038461538461, 100.0, 100.0}
  );
}
