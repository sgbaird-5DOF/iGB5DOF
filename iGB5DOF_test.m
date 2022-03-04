qm = [-0.3385 0.7952 -0.3892 -0.3187
    - 0.4534 0.5056 -0.7169 -0.1576
    - 0.4978 0.6281 -0.4914 -0.3407
    - 0.2906 0.8879 -0.3494 -0.0710
    - 0.4506 0.5342 -0.7140 0.0421];

nA = [0.2138 -0.9519 0.2195
    0.0216 -0.9924 0.1208
    0.0224 -0.9921 0.1238
    0.3663 -0.9206 -0.1353
    - 0.0681 -0.9615 0.2662];

y_brk = [1.1024
    1.0456
    1.2394
    1.0771
    1.0969];

o = [0.7756 -0.5895 -0.2066 0.0905 0.1547 0.9174 -0.3478 0.1159
    0.7289 -0.6421 -0.1652 0.1707 -0.0974 0.8081 -0.4625 0.3515
    0.7328 -0.6437 -0.1540 0.1579 0.0177 0.9108 -0.4036 0.0848
    0.6575 -0.7014 -0.2752 -0.0032 0.3353 0.8061 -0.2024 0.4436
    0.7791 -0.6003 -0.0806 0.1613 -0.0947 0.7985 -0.4085 0.4318];

omA = {
[0.512624451530477 -0.657159350320086 -0.552591856599390
    0.657159350320086 0.714497269116494 -0.240073406923707
    0.552591856599390 -0.240073406923707 0.798127182413983]

[0.121440202925350 -0.911348751314742 -0.393313781338146
    0.911348751314742 0.259384009636557 -0.319630081534770
    0.393313781338146 -0.319630081534770 0.862056193288794]

[0.371901341991074 -0.770705419971759 0.517399794599294
    0.770705419971759 0.567033848431274 0.290664360318518
    -0.517399794599294 0.290664360318518 0.804867493559800]

[0.409657635619033 -0.802217917674215 -0.434312139067595
    0.802217917674215 0.543468164768288 -0.247161417793789
    0.434312139067595 -0.247161417793789 0.866189470850745]

[0.292476840936213 -0.948627226779487 0.120680081739422
    0.948627226779487 0.303744882008488 0.088574439125010
    -0.120680081739422 0.088574439125010 0.988731958927725]
    };

omA = cat(3, omA{:});

omB = {
[0.704618012538972 -0.429422051410986 -0.564898360917793
    0.112452749791660 -0.718457626066334 0.686420438660900
    -0.700619608253953 -0.547188579293380 -0.457948494069988]

[0.808751758585647 -0.314406757225238 -0.497060342408980
    -0.041147709630021 -0.873300702443874 0.485440778162733
    -0.586709007069796 -0.372148188390252 -0.719220596827585]

[0.758865010389821 0.137595249881272 -0.636546497293142
    -0.012745190176129 -0.974100748524005 -0.225754937603474
    -0.651123226535064 0.179430429233964 -0.737457297022332]

[0.818183328375135 0.075360467241862 -0.569997229068769
    0.370203278160206 -0.827576046860801 0.421980354402532
    -0.439915416879681 -0.556272133617474 -0.705007616520729]

[0.745460588249363 -0.248582470363371 -0.618462017261334
    -0.106096989428666 -0.960278267649146 0.258087344739971
    -0.658051624267495 -0.126776985722168 -0.742223454015106]
    };

omB = cat(3, omB{:});

tic
[y_init_pred, mdl] = iGB5DOF(omA, omB, "return_model", true);
toc

save('mdl.mat', 'mdl')

tic
y_pred = iGB5DOF(omA, omB, mdl);
toc

y_pred
y_pred - y_init_pred
y_pred - y_brk
