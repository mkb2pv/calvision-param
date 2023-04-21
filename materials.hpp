// parameters are:
// 0 and 1: first and second decay times
// 2 and 3: fraction of first decay time, fraction of second decay time
TF1* p_decay_time = new TF1("p_decay","[2]/[0]*exp(-x/[0])+[3]/[1]*exp(-x/[1])",0,100);

const int nEntries_FAST = 512;
const double eV = 1.0;

double PhotonEnergy_FAST[nEntries_FAST] = {1.90769 * eV, 1.947435 * eV, 1.98718 * eV, 2.026925 * eV, 2.06667 * eV,
                                               2.09196 * eV, 2.11725 * eV, 2.14254 * eV, 2.16783 * eV, 2.16973 * eV, 2.17163 * eV, 2.17353 * eV,
                                               2.17544 * eV, 2.17735 * eV, 2.17926 * eV, 2.18118 * eV, 2.1831 * eV, 2.18502 * eV, 2.18695 * eV,
                                               2.18888 * eV, 2.19081 * eV, 2.19275 * eV, 2.19469 * eV, 2.19663 * eV, 2.19858 * eV, 2.20053 * eV,
                                               2.20249 * eV, 2.20444 * eV, 2.20641 * eV, 2.20837 * eV, 2.21034 * eV, 2.21231 * eV, 2.21429 * eV,
                                               2.21626 * eV, 2.21825 * eV, 2.22023 * eV, 2.22222 * eV, 2.22422 * eV, 2.22621 * eV, 2.22821 * eV,
                                               2.23022 * eV, 2.23222 * eV, 2.23423 * eV, 2.23625 * eV, 2.23827 * eV, 2.24029 * eV, 2.24231 * eV,
                                               2.24434 * eV, 2.24638 * eV, 2.24841 * eV, 2.25045 * eV, 2.2525 * eV, 2.25455 * eV, 2.2566 * eV,
                                               2.25865 * eV, 2.26071 * eV, 2.26277 * eV, 2.26484 * eV, 2.26691 * eV, 2.26898 * eV, 2.27106 * eV,
                                               2.27314 * eV, 2.27523 * eV, 2.27732 * eV, 2.27941 * eV, 2.28151 * eV, 2.28361 * eV, 2.28571 * eV,
                                               2.28782 * eV, 2.28994 * eV, 2.29205 * eV, 2.29417 * eV, 2.2963 * eV, 2.29842 * eV, 2.30056 * eV,
                                               2.30269 * eV, 2.30483 * eV, 2.30698 * eV, 2.30912 * eV, 2.31128 * eV, 2.31343 * eV, 2.31559 * eV,
                                               2.31776 * eV, 2.31993 * eV, 2.3221 * eV, 2.32427 * eV, 2.32645 * eV, 2.32864 * eV, 2.33083 * eV,
                                               2.33302 * eV, 2.33522 * eV, 2.33742 * eV, 2.33962 * eV, 2.34183 * eV, 2.34405 * eV, 2.34626 * eV,
                                               2.34848 * eV, 2.35071 * eV, 2.35294 * eV, 2.35518 * eV, 2.35741 * eV, 2.35966 * eV, 2.3619 * eV,
                                               2.36416 * eV, 2.36641 * eV, 2.36867 * eV, 2.37094 * eV, 2.37321 * eV, 2.37548 * eV, 2.37776 * eV,
                                               2.38004 * eV, 2.38232 * eV, 2.38462 * eV, 2.38691 * eV, 2.38921 * eV, 2.39151 * eV, 2.39382 * eV,
                                               2.39614 * eV, 2.39845 * eV, 2.40077 * eV, 2.4031 * eV, 2.40543 * eV, 2.40777 * eV, 2.41011 * eV,
                                               2.41245 * eV, 2.4148 * eV, 2.41715 * eV, 2.41951 * eV, 2.42188 * eV, 2.42424 * eV, 2.42661 * eV,
                                               2.42899 * eV, 2.43137 * eV, 2.43376 * eV, 2.43615 * eV, 2.43854 * eV, 2.44094 * eV, 2.44335 * eV,
                                               2.44576 * eV, 2.44817 * eV, 2.45059 * eV, 2.45302 * eV, 2.45545 * eV, 2.45788 * eV, 2.46032 * eV,
                                               2.46276 * eV, 2.46521 * eV, 2.46766 * eV, 2.47012 * eV, 2.47258 * eV, 2.47505 * eV, 2.47752 * eV,
                                               2.48 * eV, 2.48248 * eV, 2.48497 * eV, 2.48746 * eV, 2.48996 * eV, 2.49246 * eV, 2.49497 * eV,
                                               2.49748 * eV, 2.5 * eV, 2.50252 * eV, 2.50505 * eV, 2.50758 * eV, 2.51012 * eV, 2.51266 * eV,
                                               2.51521 * eV, 2.51777 * eV, 2.52033 * eV, 2.52289 * eV, 2.52546 * eV, 2.52803 * eV, 2.53061 * eV,
                                               2.5332 * eV, 2.53579 * eV, 2.53838 * eV, 2.54098 * eV, 2.54359 * eV, 2.5462 * eV, 2.54882 * eV,
                                               2.55144 * eV, 2.55407 * eV, 2.5567 * eV, 2.55934 * eV, 2.56198 * eV, 2.56463 * eV, 2.56729 * eV,
                                               2.56995 * eV, 2.57261 * eV, 2.57529 * eV, 2.57796 * eV, 2.58065 * eV, 2.58333 * eV, 2.58603 * eV,
                                               2.58873 * eV, 2.59143 * eV, 2.59414 * eV, 2.59686 * eV, 2.59958 * eV, 2.60231 * eV, 2.60504 * eV,
                                               2.60778 * eV, 2.61053 * eV, 2.61328 * eV, 2.61603 * eV, 2.6188 * eV, 2.62156 * eV, 2.62434 * eV,
                                               2.62712 * eV, 2.6299 * eV, 2.6327 * eV, 2.63549 * eV, 2.6383 * eV, 2.64111 * eV, 2.64392 * eV,
                                               2.64674 * eV, 2.64957 * eV, 2.65241 * eV, 2.65525 * eV, 2.65809 * eV, 2.66094 * eV, 2.6638 * eV,
                                               2.66667 * eV, 2.66954 * eV, 2.67241 * eV, 2.6753 * eV, 2.67819 * eV, 2.68108 * eV, 2.68398 * eV,
                                               2.68689 * eV, 2.6898 * eV, 2.69273 * eV, 2.69565 * eV, 2.69859 * eV, 2.70153 * eV, 2.70447 * eV,
                                               2.70742 * eV, 2.71038 * eV, 2.71335 * eV, 2.71632 * eV, 2.7193 * eV, 2.72228 * eV, 2.72527 * eV,
                                               2.72827 * eV, 2.73128 * eV, 2.73429 * eV, 2.73731 * eV, 2.74033 * eV, 2.74336 * eV, 2.7464 * eV,
                                               2.74945 * eV, 2.7525 * eV, 2.75556 * eV, 2.75862 * eV, 2.76169 * eV, 2.76477 * eV, 2.76786 * eV,
                                               2.77095 * eV, 2.77405 * eV, 2.77716 * eV, 2.78027 * eV, 2.78339 * eV, 2.78652 * eV, 2.78965 * eV,
                                               2.79279 * eV, 2.79594 * eV, 2.7991 * eV, 2.80226 * eV, 2.80543 * eV, 2.80861 * eV, 2.81179 * eV,
                                               2.81498 * eV, 2.81818 * eV, 2.82139 * eV, 2.8246 * eV, 2.82782 * eV, 2.83105 * eV, 2.83429 * eV,
                                               2.83753 * eV, 2.84078 * eV, 2.84404 * eV, 2.8473 * eV, 2.85057 * eV, 2.85386 * eV, 2.85714 * eV,
                                               2.86044 * eV, 2.86374 * eV, 2.86705 * eV, 2.87037 * eV, 2.8737 * eV, 2.87703 * eV, 2.88037 * eV,
                                               2.88372 * eV, 2.88708 * eV, 2.89044 * eV, 2.89382 * eV, 2.8972 * eV, 2.90058 * eV, 2.90398 * eV,
                                               2.90739 * eV, 2.9108 * eV, 2.91422 * eV, 2.91765 * eV, 2.92108 * eV, 2.92453 * eV, 2.92798 * eV,
                                               2.93144 * eV, 2.93491 * eV, 2.93839 * eV, 2.94187 * eV, 2.94537 * eV, 2.94887 * eV, 2.95238 * eV,
                                               2.9559 * eV, 2.95943 * eV, 2.96296 * eV, 2.96651 * eV, 2.97006 * eV, 2.97362 * eV, 2.97719 * eV,
                                               2.98077 * eV, 2.98436 * eV, 2.98795 * eV, 2.99156 * eV, 2.99517 * eV, 2.99879 * eV, 3.00242 * eV,
                                               3.00606 * eV, 3.00971 * eV, 3.01337 * eV, 3.01703 * eV, 3.02071 * eV, 3.02439 * eV, 3.02808 * eV,
                                               3.03178 * eV, 3.0355 * eV, 3.03922 * eV, 3.04294 * eV, 3.04668 * eV, 3.05043 * eV, 3.05419 * eV,
                                               3.05795 * eV, 3.06173 * eV, 3.06551 * eV, 3.06931 * eV, 3.07311 * eV, 3.07692 * eV, 3.08075 * eV,
                                               3.08458 * eV, 3.08842 * eV, 3.09227 * eV, 3.09613 * eV, 3.1 * eV, 3.10388 * eV, 3.10777 * eV, 3.11167 * eV,
                                               3.11558 * eV, 3.1195 * eV, 3.12343 * eV, 3.12736 * eV, 3.13131 * eV, 3.13527 * eV, 3.13924 * eV, 3.14322 * eV,
                                               3.14721 * eV, 3.15121 * eV, 3.15522 * eV, 3.15924 * eV, 3.16327 * eV, 3.16731 * eV, 3.17136 * eV, 3.17542 * eV,
                                               3.17949 * eV, 3.18357 * eV, 3.18766 * eV, 3.19176 * eV, 3.19588 * eV, 3.2 * eV, 3.20413 * eV, 3.20828 * eV,
                                               3.21244 * eV, 3.2166 * eV, 3.22078 * eV, 3.22497 * eV, 3.22917 * eV, 3.23338 * eV, 3.2376 * eV, 3.24183 * eV,
                                               3.24607 * eV, 3.25033 * eV, 3.25459 * eV, 3.25887 * eV, 3.26316 * eV, 3.26746 * eV, 3.27177 * eV, 3.27609 * eV,
                                               3.28042 * eV, 3.28477 * eV, 3.28912 * eV, 3.29349 * eV, 3.29787 * eV, 3.30226 * eV, 3.30667 * eV, 3.31108 * eV,
                                               3.31551 * eV, 3.31995 * eV, 3.3244 * eV, 3.32886 * eV, 3.33333 * eV, 3.33782 * eV, 3.34232 * eV, 3.34683 * eV,
                                               3.35135 * eV, 3.35589 * eV, 3.36043 * eV, 3.36499 * eV, 3.36957 * eV, 3.37415 * eV, 3.37875 * eV, 3.38336 * eV,
                                               3.38798 * eV, 3.39261 * eV, 3.39726 * eV, 3.40192 * eV, 3.40659 * eV, 3.41128 * eV, 3.41598 * eV, 3.42069 * eV,
                                               3.42541 * eV, 3.43015 * eV, 3.4349 * eV, 3.43967 * eV, 3.44444 * eV, 3.44924 * eV, 3.45404 * eV, 3.45886 * eV,
                                               3.46369 * eV, 3.46853 * eV, 3.47339 * eV, 3.47826 * eV, 3.48315 * eV, 3.48805 * eV, 3.49296 * eV, 3.49788 * eV,
                                               3.50282 * eV, 3.50778 * eV, 3.51275 * eV, 3.51773 * eV, 3.52273 * eV, 3.52774 * eV, 3.53276 * eV, 3.5378 * eV,
                                               3.54286 * eV, 3.54793 * eV, 3.55301 * eV, 3.55811 * eV, 3.56322 * eV, 3.56835 * eV, 3.57349 * eV, 3.57864 * eV,
                                               3.58382 * eV, 3.589 * eV, 3.5942 * eV, 3.59942 * eV, 3.60465 * eV, 3.6099 * eV, 3.61516 * eV, 3.62044 * eV,
                                               3.62573 * eV, 3.63104 * eV, 3.63636 * eV, 3.6417 * eV, 3.64706 * eV, 3.65243 * eV, 3.65782 * eV, 3.66322 * eV,
                                               3.66864 * eV, 3.67407 * eV, 3.67953 * eV, 3.68499 * eV, 3.69048 * eV, 3.69598 * eV, 3.70149 * eV, 3.70703 * eV,
                                               3.71257 * eV, 3.71814 * eV, 3.72372 * eV, 3.72932 * eV, 3.73494 * eV, 3.74057 * eV, 3.74622 * eV, 3.75189 * eV,
                                               3.75758 * eV, 3.76328 * eV, 3.769 * eV, 3.77473 * eV, 3.78049 * eV, 3.78626 * eV, 3.79205 * eV, 3.79786 * eV,
                                               3.80368 * eV, 3.80952 * eV, 3.81538 * eV, 3.82126 * eV, 3.82716 * eV, 3.83308 * eV, 3.83901 * eV, 3.84496 * eV,
					   3.85093 * eV, 3.85692 * eV, 3.86293 * eV, 3.86895 * eV};




double FastComponent[nEntries_FAST] = {0.00052381, 0.00420238, 0.00788096,
                                           0.01155953, 0.0152381, 0.01877668, 0.02231525, 0.02585382, 0.0293924, 0.0269467,
                                           0.0256848, 0.0254981, 0.0260076, 0.0268133, 0.0273648, 0.0280457, 0.0287067, 0.0289924,
                                           0.029, 0.0288914, 0.0289476, 0.0291505, 0.0292914, 0.0290752, 0.0287514, 0.0286714,
                                           0.0288343, 0.0293133, 0.0302143, 0.031561, 0.0330762, 0.034459, 0.0351705, 0.03478,
                                           0.0340324, 0.0329552, 0.0314514, 0.0300962, 0.0292152, 0.028841, 0.0289524, 0.0294419,
                                           0.0297876, 0.0303133, 0.0310962, 0.0317743, 0.0322238, 0.0325276, 0.0329733, 0.03336,
                                           0.0337895, 0.0342257, 0.0347362, 0.0354695, 0.0363, 0.0369238, 0.036979, 0.036579,
                                           0.0356076, 0.0343095, 0.032999, 0.0316876, 0.0306171, 0.0302533, 0.0308752, 0.0322067,
                                           0.0339943, 0.0360067, 0.0384286, 0.0411333, 0.04334, 0.0450924, 0.0462695, 0.0472533,
                                           0.0481019, 0.0483352, 0.0482181, 0.04864, 0.0499019, 0.0517543, 0.0543505, 0.0575267,
                                           0.0607876, 0.0641314, 0.0667838, 0.0683514, 0.0693419, 0.0702543, 0.0710981, 0.0720552,
                                           0.0736676, 0.0752762, 0.0773286, 0.0791752, 0.0807333, 0.082079, 0.0833629, 0.0845933,
                                           0.0859524, 0.0877581, 0.0892943, 0.0910914, 0.0929019, 0.0952905, 0.0978371, 0.100682,
                                           0.103224, 0.105718, 0.107852, 0.109241, 0.109795, 0.109354, 0.109503, 0.110292, 0.112409,
                                           0.115256, 0.11961, 0.124909, 0.130732, 0.13613, 0.140314, 0.144119, 0.14665, 0.148469,
                                           0.149407, 0.150265, 0.151147, 0.152253, 0.154354, 0.15737, 0.161138, 0.164838, 0.168873,
                                           0.172971, 0.177416, 0.181659, 0.185039, 0.18823, 0.191871, 0.195512, 0.198159, 0.20033,
                                           0.202124, 0.203936, 0.205858, 0.207766, 0.20994, 0.212648, 0.216484, 0.22057, 0.22497,
                                           0.229337, 0.232884, 0.23596, 0.239326, 0.242921, 0.245595, 0.248666, 0.2521, 0.255899,
                                           0.260253, 0.264172, 0.268089, 0.27221, 0.27688, 0.280961, 0.284639, 0.28845, 0.291873,
                                           0.294659, 0.29668, 0.298346, 0.300048, 0.302431, 0.305583, 0.309162, 0.313943, 0.319651,
                                           0.324986, 0.329565, 0.333289, 0.336565, 0.33924, 0.341681, 0.343442, 0.345074, 0.347293,
                                           0.349824, 0.352279, 0.354626, 0.357285, 0.360091, 0.363057, 0.366477, 0.37011, 0.373876,
                                           0.377648, 0.380938, 0.384106, 0.387288, 0.390137, 0.39199, 0.393832, 0.39634, 0.399255,
                                           0.402285, 0.405265, 0.408943, 0.412963, 0.417683, 0.42203, 0.425902, 0.430032, 0.434444,
                                           0.438648, 0.442674, 0.445648, 0.44771, 0.450105, 0.452818, 0.455475, 0.459257, 0.46497,
                                           0.471776, 0.4804, 0.489469, 0.497992, 0.507131, 0.516506, 0.525018, 0.532744, 0.540899,
                                           0.548129, 0.555484, 0.562471, 0.568386, 0.574721, 0.581799, 0.58892, 0.595488, 0.603247,
                                           0.611145, 0.619414, 0.628492, 0.637634, 0.646924, 0.656813, 0.66681, 0.675362, 0.683059,
                                           0.689791, 0.69601, 0.702039, 0.708148, 0.714768, 0.722988, 0.732344, 0.74167, 0.75073,
                                           0.759785, 0.768961, 0.777719, 0.78556, 0.792373, 0.79914, 0.80549, 0.81113, 0.815989,
                                           0.821, 0.826681, 0.832696, 0.838765, 0.844266, 0.850464, 0.85645, 0.862316, 0.86791,
                                           0.873631, 0.879935, 0.887078, 0.895343, 0.903276, 0.910896, 0.917349, 0.922744, 0.927535,
                                           0.931457, 0.93481, 0.938199, 0.942202, 0.946601, 0.950905, 0.954545, 0.957467, 0.96089,
                                           0.964401, 0.967275, 0.968981, 0.97016, 0.971976, 0.973287, 0.973908, 0.973841, 0.97475,
                                           0.977297, 0.981698, 0.986047, 0.990052, 0.994456, 0.998224, 1.0007, 1.00106, 1.00068,
                                           1.00017, 1.00024, 1.00012, 0.99895, 0.997525, 0.995587, 0.994011, 0.991989, 0.990154,
                                           0.988663, 0.987812, 0.988681, 0.990095, 0.990958, 0.990571, 0.990366, 0.989786, 0.988225,
                                           0.985411, 0.981308, 0.977582, 0.973715, 0.96889, 0.963269, 0.958067, 0.954727, 0.952326,
                                           0.95059, 0.949261, 0.949456, 0.949965, 0.948949, 0.946274, 0.94231, 0.937687, 0.93214,
                                           0.926558, 0.920441, 0.915683, 0.912037, 0.908574, 0.904968, 0.901723, 0.898805, 0.895016,
                                           0.891226, 0.886231, 0.880648, 0.874508, 0.867607, 0.859656, 0.851617, 0.844196, 0.83623,
                                           0.828706, 0.822149, 0.817083, 0.812778, 0.808321, 0.803222, 0.798333, 0.793735, 0.787804,
                                           0.780487, 0.772463, 0.764901, 0.75783, 0.750741, 0.743238, 0.737445, 0.73221, 0.725617,
                                           0.717075, 0.707011, 0.696076, 0.684175, 0.670404, 0.65462, 0.640174, 0.627405, 0.6152,
                                           0.603558, 0.592237, 0.58155, 0.57139, 0.559997, 0.546431, 0.532181, 0.517833, 0.503294,
                                           0.488553, 0.474083, 0.460749, 0.449591, 0.439908, 0.431058, 0.42282, 0.414699, 0.406633,
                                           0.398634, 0.39069, 0.382162, 0.373201, 0.364355, 0.355435, 0.346777, 0.337376, 0.32759,
                                           0.31762, 0.307741, 0.29736, 0.286301, 0.274514, 0.262031, 0.249232, 0.235434, 0.220187,
                                           0.204732, 0.18955, 0.174084, 0.158599, 0.143255, 0.128051, 0.113468, 0.0996657, 0.0862448,
                                           0.0741762, 0.0637238, 0.054501, 0.0470933, 0.0412562, 0.0365495, 0.0324981, 0.0291943,
                                           0.0259467, 0.0229876, 0.0201476, 0.0172495, 0.0144133, 0.0121181, 0.010861, 0.0100343,
                                           0.00974476, 0.0103733, 0.0119886, 0.01364, 0.0151286, 0.0161257, 0.0168276, 0.0176267,
                                           0.0177667, 0.0169867, 0.01598, 0.015241, 0.0144143, 0.0135886, 0.0125457, 0.0115524,
                                           0.0113305, 0.0114295, 0.0114038, 0.0114352, 0.01208, 0.0132114, 0.0141905, 0.0147667,
                                           0.0149648, 0.0148695, 0.0140505, 0.0127952, 0.0109514, 0.00864667, 0.00670762, 0.00527143,
                                           0.0046019, 0.00473524, 0.00552476, 0.0065, 0.00768667, 0.0084381, 0.00831333, 0.00752286,
                                           0.0062181, 0.00454952, 0.00287905, 0.00136476, 0.000487619, 0.000487619, 0.000514286,
                                           0.000467619, 0.000337143, 0.00047619, 0.00104, 0.00124, 0.000652381, 0.0015, 0.00581905,
				       0.0120495, 0.0200286};


const int nEntries_SCY = 12;
const double MeV = 1;

double ElectronEnergy_SCY[nEntries_SCY] = {0.000 * MeV, 0.015 * MeV, 0.020 * MeV, 0.030 * MeV,
       0.040 * MeV, 0.060 * MeV, 0.080 * MeV, 0.090 * MeV,
					   0.105 * MeV, 0.300 * MeV, 0.500 * MeV, 1.000 * MeV};

double ScintilYield[nEntries_SCY] = {0.10, 0.46, 0.60, 0.68,
       0.74, 0.80, 0.82, 0.84,
				     0.87, 0.96, 0.98, 1.00};
