function varargout = iGB5DOF(P, Q, mdl, nv)

    arguments
        P(3, 3, :) double
        Q(3, 3, :) double
        mdl = []
        nv.return_model(1, 1) logical = false
        nv.numNearest(1, 1) {mustBeInteger} = 50
    end

    %  GB5DOF  computes the energy of an arbitrary boundary in FCC metals (BRK energy function)
    %
    %     en = GB5DOF(P,Q) computes the energy of a boundary
    %     described by two rotation matrices P and Q. The character string
    %     variable AlCuParameter can take on one of four values: 'Al', 'Ni',
    %     'Au', or 'Cu'.  The output variable en is the computed boundary energy
    %     in units J/m^2.
    %
    %     P and Q are two properly normalized 3x3 rotation matrices defining
    %     orientations of the two grains with respect to a laboratory (sample)
    %     frame. For any vector V expressed in the cube frame of grain P (or Q),
    %     P*V (or Q*V) expresses the same vector in the laboratory frame. By
    %     convention, the first row of P (or Q) is the row vector of the boundary
    %     plane normal Np = P(1,:) (or Nq = Q(1,:)) written in the cube frame
    %     of grain P (or Q). Thus, P*Np' = Q*Nq' = [1 0 0]'.
    %
    %     Examples
    %
    %       With P and Q matrices defined as follows
    %
    %           P = [ 0.5774    0.5774    0.5774 ;
    %                 0.7071   -0.7071         0 ;
    %                 0.4082    0.4082   -0.8165 ]
    %
    %           Q = [ 0.5774    0.5774    0.5774 ;
    %               [-0.7071    0.7071         0 ;
    %                -0.4082   -0.4082    0.8165 ]
    return_model = nv.return_model;

    if isempty(mdl)
        [qm, nA, y] = get_olmsted_ni_data();
        numNearest = nv.numNearest;
        [~, ~, mdl, ~] = interp5DOF(qm, nA, y, 'nNN', numNearest, 'mygpropts', {'PredictMethod', 'exact', 'Sigma', 0.01});
    end

    npts = size(Q, 3);
    o2 = zeros(npts, 8);

    for i = 1:npts
        o2(i, :) = om2oct(P(:, :, i), Q(:, :, i));
    end

    [qm2, nA2] = oct2five(o2);
    y_pred = mdl.interpfn(qm2, nA2);

    varargout{1} = y_pred;

    if return_model
        varargout{2} = mdl;
    end

end

function [ypred, interpfn, mdl, mdlpars] = interp5DOF(qm, nA, y, qm2, nA2, method, epsijk, nv)

    arguments
        qm %input misorientation quaternions
        nA %input BP normals
        y(:, 1) %property values
        qm2 = [] %query misorientations
        nA2 = [] %query BP normals
        method char {mustBeMember(method, {'gpr', 'sphgpr', 'pbary', 'sphbary', 'idw', 'nn', 'avg'})} = 'gpr'
        epsijk(1, 1) double = 1
        nv.pgnum(1, 1) double = 32 %m-3m (i.e. m\overbar{3}m) FCC symmetry default
        nv.databary = [] %for use with bary methods
        nv.facetIDs = [] %for use with bary methods
        nv.ytrue = [] %user-specified "true" values for error calculations
        nv.modelparsspec = struct()
        nv.brkQ(1, 1) logical = false %whether to compute BRK values as ytrue
        nv.sigma(1, 1) double = 0 %noise to add to property values
        nv.mygpropts = struct.empty %for use with gpr methods 'gpr' or 'sphgpr'
        nv.r double = [] %for use with 'idw' method, alternatively set to [] for automatic estimation
        nv.uuid(1, 8) char = get_uuid() %unique ID associated with this interpolation run
        nv.o = [] %input octonions, specify these or qm/nA pairs
        nv.o2 = [] %query octonions, specify these or qm2/nA2 pairs
        nv.oref = get_ocubo(1, 'random', [], 10)
        nv.nforce double = 1
        nv.nforceQ(1, 1) logical = false
        nv.dispQ(1, 1) logical = true
        nv.IncludeTies(1, 1) {mustBeLogical} = true
        nv.nNN(1, 1) double = 1
        nv.projQ(1, 1) logical = true
    end

    % INTERP5DOF  Convert misorientation and boundary plane normal 5DOF input
    % data to a closed, octonion, hyperspherical mesh and interpolate property
    % values for arbitrary grain boundaries using spherical barycentric
    % coordinates, planar barycentric coordinates, or a Gaussian process
    % regression.
    %--------------------------------------------------------------------------
    % Inputs:
    %  qm - list of misorientation quaternions (data), as in qmult(qB,qinv(qA))
    %  for grains A and B, where qA and qB are in the sample frame (same for
    %  qm2)
    %
    %  nA - list of boundary plane Cartesian unit normals (grain A frame)
    %
    %  propList - property value for each grain boundary (GB)
    %
    %  qm2 - list of misorientation quaternions for query GBs
    %
    %  nA2 - list of boundary plane Cartesian unit normals for query GBs (grain
    %   A frame)
    %
    %  method - interpolation scheme to use. Possible methods are:
    %    'sphbary' - spherical barycentric coordinates
    %    'pbary' - planar barycentric coordinates
    %    'gpr' - Gaussian process regression
    %    'nn' - nearest neighbor interpolation
    %    'insertnamehere' - no functionality, but contains instructions for
    %      implementing your own interpolation scheme
    %
    %  nv - method-specific name-value pairs
    %     method == 'sphbary' or 'pbary'
    %       'databary' - supply barycentric coordinates to reduce redundant
    %       computation if repeating interp5DOF calls for same list of qm/nA
    %       pairs, even if property values are different. I.e. the interpolation
    %       is calculated via a simple table lookup and a dot product.
    %       facetprops must also be supplied.
    %
    %       'facetprops' - supply properties of each facet vertex (i.e. same
    %       size as databary) for barycentric data interpolation. databary must
    %       also be supplied
    %
    %       'brkQ' - logical, whether or not to calculate BRK energy values for
    %       the query points to use for error calculations. If false and
    %       ytrue is not supplied, then data.props is assigned a NaN vector
    %       of the same size as mesh.props
    %
    %       'ytrue' - user supplied properties for query points for error
    %       calculations. If not supplied, data.props depends on brkQ
    %
    %       'modelparsspec' - user supplied struct of model-specific
    %       parameters, only gets used if databary and facetprops also supplied
    %       via name-value pairs
    %
    %     method == 'gpr'
    %       'gpropts' - options structure that will be passed to fitrgp() in place
    %       of the defaults that are supplied later in this function.
    %
    % Outputs:
    %  propOut - interpolated property values of queried grain boundaries
    %
    %  varargout
    %   method == 'gpr'
    %    varargout{1} - gprMdl, a Guassian Process Regression Object
    %    varargout{2} - ysd, standard deviations of predicted values (propOut)
    %    varargout{3} - yint, 95% confidence intervals of predicted values
    %
    %   method == 'sphbary' or 'pbary'
    %    varargout{1} - databary, barycentric coordinates
    %    varargout{2} - fname, which contains get_interp() workspace
    %
    % Usage:
    %  propOut = interp5DOF(qm,nA,propList,qm2,nA2,'gpr')
    %  propOut = interp5DOF(qm,nA,propList,qm2,nA2,'sphbary')
    %  propOut = interp5DOF(qm,nA,propList,qm2,nA2,'pbary')
    %
    %  [propOut,gprMdl] = interp5DOF(qm,nA,propList,qm2,nA2,'gpr')
    %  [propOut,gprMdl,ysd,yint] = interp5DOF(qm,nA,propList,qm2,nA2,'gpr')
    %
    %  [propOut] = interp5DOF(qm,nA,propList,qm2,nA2,'gpr','gpropts',gpropts);
    %
    %  [propOut,databary,fname] = interp5DOF(qm,nA,propList,qm2,nA2,'sphbary')
    %  [propOut,databary,fname] = interp5DOF(qm,nA,propList,qm2,nA2,'pbary')
    %
    % Simple Example Data
    %  npts = 100;
    %  qm = get_cubo(npts); nA = normr(rand(npts,3)); %random (qm,nA) pairs
    %  propList = rand(npts); %random property values
    %  qm2 = get_cubo(npts); nA2 = normr(rand(npts,3)); %random (qm,nA) pairs
    %  %(I suggest you look at interp5DOF_test.m instead)
    %
    % Dependencies:
    %  MATLAB 2019b or higher (mainly for the "arguments" syntax checking at
    %  the beginning of functions, which is used extensively throughout)
    %
    %  Toolboxes
    %  -Statistics and Machine Learning Toolbox
    %  -Symbolic Math Toolbox (optional, for numStabBary.m)
    %
    %  Functions
    %   see "FUNCTION DEPENDENCIES" section at end of this file (prior to "CODE
    %   GRAVEYARD" section)
    %
    % Notes:
    %  Simpler, plug & play, input/output version of run.m (different options)
    %
    %  Dependencies updated 2020-09-03 SGB
    %
    %  *If you have the Computer Vision toolbox, normr.m will be shadowed by
    %  the corresponding function in the toolbox. However, both should produce
    %  the same results for the purposes here.
    %
    %  **If you have addpathdir.m available on your path, all the other
    %  dependencies should be added as long as the functions are in a
    %  sub-folder of your current working directory. Alternatively, call the
    %  following line of code while in the "code" directory, and then move
    %  to your directory of choice:
    %  addpathdir({'normr.m','GB5DOF_setup.m','cu2qu.m','q2rod.m','five2oct.m','correctdis.m'})
    %
    %  In the context of this function, mesh is equivalent to predictors &
    %  predictor responses, and data is equivalent to the query points you want
    %  to interpolate or predict at.
    %
    %  Test functions (e.g. five2oct_test.m, get_octpairs_test.m) are
    %  available for many of the functions listed in "FUNCTION DEPENDENCIES" at
    %  the end which can help with understanding and visualizing what each
    %  function does. These are also very useful for debugging. The test
    %  functions may have other dependencies than the ones listed in this file,
    %  such as plot5DOF.m. There is also a "test" function for this function
    %  (interp5DOF.m), namely interp5DOF_test.m
    %
    %  If you want to implement a custom interpolation scheme other than the
    %  three available here, see the instructions beneath case 'insertnamehere'
    %  in the switch method statement.
    %
    %  You can minimize this preamble text in MATLAB by clicking the "minus"
    %  symbol at the top-left.
    %
    % Author(s): Sterling Baird
    %
    % Date: 2020-09-03
    %--------------------------------------------------------------------------

    %unpack (some) name-value pairs
    pgnum = nv.pgnum;
    brkQ = nv.brkQ;
    uuid = nv.uuid;
    ytrue = nv.ytrue;
    sigma = nv.sigma;
    nforceQ = nv.nforceQ;
    nforce = nv.nforce;
    dispQ = nv.dispQ;
    o2 = nv.o2;
    IncludeTies = nv.IncludeTies;
    nNN = nv.nNN;
    projQ = nv.projQ;

    %display method
    if dispQ
        disp(['method = ' method])
    end

    % add relevant folders to path (by searching subfolders for functions)
    addpath(genpath('.'))
    % addpathdir({'normr.m','GB5DOF_setup.m','cu2qu.m','q2rod.m','five2oct.m',...
    %     'correctdis.m','interp_gpr.m'})

    %% convert to octonions & symmetrize
    tic
    %predictor points
    if isempty(qm) && isempty(nA) && ~isempty(nv.o)
        predinput = 'octonion';
        otmp = nv.o;
    else
        predinput = '5dof';
        otmp = five2oct(qm, nA, epsijk);
    end

    %query points
    if isempty(qm2) && isempty(nA2) && ~isempty(o2)
        queryinput = 'octonion';
        otmp2 = o2;
    elseif ~isempty(qm2) && ~isempty(nA2) && isempty(o2)
        queryinput = '5dof';
        otmp2 = five2oct(qm2, nA2, epsijk);
    elseif isempty(qm2) && isempty(nA2) && isempty(o2)
        queryinput = 'octonion';
        otmp2 = get_ocubo(1, 'random', [], []); %dummy variable
    end

    %symmetrization
    % wtol = 1e-6;
    [o, oref, ~, ids] = get_octpairs(otmp, 'pgnum', pgnum, 'oref', nv.oref, 'IncludeTies', IncludeTies, 'nNN', nNN);
    ninputpts = size(o, 1);

    %symmetrization
    [o2, oref2] = get_octpairs(otmp2, 'pgnum', pgnum, 'oref', nv.oref, 'IncludeTies', false, 'nNN', 1);
    npredpts = size(o2, 1);

    %make sure that reference octonions are identical within tolerance
    if ~ismembertol(oref, oref2, 'ByRows', true)
        disp(['oref  == ' num2str(oref)])
        disp(['oref2 == ' num2str(oref2)])
        warning('oref ~= oref2')
    end

    symruntime = toc;

    [~, ~, nnmu, nnsigma] = get_knn(o, 'omega', 1);

    if dispQ
        disp(['ninputpts = ' int2str(ninputpts) ', npredpts = ' int2str(npredpts)])
    end

    %% projection
    %important to project both sets together so they use same SVD decomposition

    projtol = 1e-4;
    zeroQ = false;
    o = normr(o);
    o2 = normr(o2);
    [a, usv] = proj_down([o; o2], projtol, 'zero', zeroQ, 'force', nforceQ, 'nforcedim', nforce);

    d = size(a, 2);
    %projected points
    if d <= 8
        ppts = proj_down(o, projtol, usv, 'zero', zeroQ, 'force', nforceQ, 'nforcedim', nforce);
        ppts2 = proj_down(o2, projtol, usv, 'zero', zeroQ, 'force', nforceQ, 'nforcedim', nforce);
    else
        error("Input doesn't have degenerate dimension or has too few (i.e. check input data), or try reducing proj_down tolerance input (tol)")
    end

    %% mesh and data struct setup
    %mesh
    mesh.pts = o;
    mesh.ppts = ppts;
    mesh.npts = ninputpts;

    %data
    data.pts = o2;
    data.ppts = ppts2;
    data.npts = npredpts;

    %mesh property values
    if brkQ

        if ~isempty(y)
            warning('overriding "y" values with BRK values')
        end

        pA = o(:, 1:4);
        pB = o(:, 5:8);
        mA = [0 0 1];
        y = GB5DOF_setup(pA, pB, mA, 'Ni', epsijk);

        % add noise to BRK values
        noisetype = 'normal';

        switch noisetype
            case 'normal'
                %             y = normrnd(y,sigma);
                y = y + abs(normrnd(zeros(size(y)), sigma));
            case 'uniform'
                y = y + sigma * 2 * (rand(size(y)) - 0.5); %uniform
        end

    elseif isempty(y)
        y = nan(size(ppts2, 1), 1);
    end

    y = y(ids); %data augmentation for property values

    mesh.props = y;

    %data property values
    if brkQ

        if ~isempty(ytrue)
            warning('overriding "ytrue" values with BRK values')
        end

        pA = o2(:, 1:4);
        pB = o2(:, 5:8);
        mA = [0 0 1];
        ytrue = GB5DOF_setup(pA, pB, mA, 'Ni', epsijk);

        % add noise to BRK values if applicable
        noisetype = 'normal';

        switch noisetype
            case 'normal'
                %             ytrue = normrnd(ytrue,sigma);
                ytrue = ytrue + abs(normrnd(zeros(size(ytrue)), sigma));
            case 'uniform'
                ytrue = ytrue + sigma * 2 * (rand(size(ytrue)) - 0.5); %uniform
        end

    elseif isempty(ytrue)
        ytrue = nan(size(ppts2, 1), 1);
    end

    data.props = ytrue;

    %% additional variables
    % current date and time
    starttime = datetime(clock);
    % number of cores (i.e. parfor workers)
    p = gcp;
    ncores = p.NumWorkers;
    %git commit version
    gitcommit = get_gitcommit();

    %% package into struct
    %general model variables
    mdlgen = var_names(method, projtol, zeroQ, usv, starttime, ncores, ninputpts, npredpts, ...
    gitcommit, uuid, predinput, queryinput, projQ, oref, oref2, nnmu, nnsigma, symruntime, ...
        IncludeTies, nNN);
    %general parameters
    mdlparsgen = var_names(method, projtol, zeroQ, starttime, ninputpts, ...
    npredpts, ncores, gitcommit, uuid, predinput, queryinput, projQ, oref, oref2, nnmu, nnsigma, ...
        symruntime, IncludeTies, nNN);

    %% method-specific interpolation
    tic

    switch method
        case {'sphbary', 'pbary'}

            if isempty(nv.databary)
                %% get triangulation
                [pptstmp, usvtri] = proj_down(o, projtol, 'zero', true);
                K = sphconvhulln(pptstmp);
                mesh.pts = proj_up(pptstmp, usvtri);

                %% compute intersecting facet IDs
                nnMax = 10;
                inttol = 0.01;
                disp('intersect_facet')
                intfacetIDs = intersect_facet(ppts, K, ppts2, inttol, 'inttype', 'planar', 'nnMax', nnMax);

                %% mesh triangulation and filename
                mesh.sphK = K;

                %% interpolation
                disp('interpolation')
                %method-specific parameters
                switch method
                    case 'sphbary'
                        barytol = 0.2;
                        barytype = 'spherical';
                        mesh.ppts = normr(mesh.ppts);
                        data.ppts = normr(data.ppts);

                    case 'pbary'
                        barytol = 0.1;
                        barytype = 'planar';
                end

                %interpolation
                [ypred, databary, facetprops, facetIDs, barypars] = get_interp(mesh, data, intfacetIDs, barytype, barytol);

                %model command and interpolation function
                mdlcmd = @(mesh, data, intfacetIDs, barytype, barytol) get_interp(mesh, data, intfacetIDs, barytype, barytol);
                interpfn = @(qm2, nA2) interp_bary(mesh, [], qm2, nA2, usv, zeroQ, barytype, barytol, projtol, nnMax, brkQ);

                %unpack intersection metrics
                nints = barypars.nints;
                numnonints = barypars.numnonints;
                int_fraction = barypars.int_fraction;

                %unpack NN extrapolation RMSE and MAE values
                nn_rmse = barypars.nn_errmetrics.rmse;
                nn_mae = barypars.nn_errmetrics.mae;

                %model-specific variables
                mdlspec = var_names(databary, facetprops, barytol, barytype, inttol, ...
                intfacetIDs, nnMax, facetIDs, barypars, nn_rmse, nn_mae);

                %model-specific parameters
                mdlparsspec = var_names(barytol, barytype, inttol, nnMax, ...
                nn_rmse, nn_mae, nints, numnonints, int_fraction, barypars);

            else
                %unpack
                databary = nv.databary;
                facetIDs = nv.facetIDs;

                %interpolate using supplied barycentric coordinates
                [ypred, facetprops, NNextrapID, nnList] = interp_bary_fast(ppts, ppts2, meshprops, databary, facetIDs);

                %model command and interp function
                mdlcmd = @(databary, facetprops) dot(databary, facetprops, 2);
                interpfn = @(propList) interp_bary_fast(ppts, ppts2, meshprops, databary, facetIDs);

                %model-specific variables
                mdlspec = var_names(databary, facetprops, NNextrapID, nnList, facetIDs);

                %model-specific parameters
                mdlparsspec = nv.modelparsspec;
            end

        case {'sphgpr', 'gpr'}

            if projQ
                X = ppts;
                X2 = ppts2;
            else
                X = o;
                X2 = o2;
            end

            %gpr options
            if isempty(nv.mygpropts)
                %% interp5DOF's default gpr options
                thresh = Inf;

                if ninputpts <= thresh
                    PredictMethod = 'fic';
                else
                    PredictMethod = 'bcd';
                    gpropts = {'BlockSize', 10000};
                end

                if exist('gpropts', 'var') ~= 1
                    gpropts = {};
                end

                gpropts = [gpropts {'PredictMethod', PredictMethod}];
                %             gpropts = {'PredictMethod',PredictMethod};

                if strcmp(method, 'sphgpr')
                    %squared exponential kernel function with octonion distance
                    kfcn = @(XN, XM, theta) (exp(theta(2))^2) * exp(- (pdist2(XN, XM, @get_alen).^2) / (2 * exp(theta(1))^2));
                    theta0 = [deg2rad(10), std(y) / sqrt(2)]; %initial length scale and noise
                    gpropts = [gpropts, {'KernelFunction', kfcn, 'KernelParameters', theta0}];
                end

            else
                % user-supplied gpr options
                gpropts = nv.mygpropts;
                numopts = length(gpropts);
                gproptnames = gpropts(1:2:numopts);
                gproptvals = gpropts(2:2:end);
                gproptstruct = cell2struct(gproptvals, gproptnames, 2);

                %extract parameters (for table)
                G = gproptstruct;

                if isempty(fieldnames(G))
                    error('user-specified gpropts is empty')
                end

                gproptshort = struct();

                if isfield(G, 'HyperparameterOptimizationOptions')
                    G1 = G.HyperparameterOptimizationOptions;

                    if isfield(G1, 'UseParallel')
                        gprParallelQ = G1.UseParallel;
                        gproptshort.gprParallelQ = gprParallelQ;
                    end

                    if isfield(G1, 'Optimizer')
                        hyperoptimizer = G1.Optimizer;
                        gproptshort.hyperoptimizer = hyperoptimizer;
                    end

                    if isfield(G1, 'MaxObjectiveEvaluations')
                        maxhyperobj = G1.MaxObjectiveEvaluations;
                        gproptshort.maxhyperobj = maxhyperobj;
                    end

                end

                if isfield(G, 'PredictMethod')
                    PredictMethod = G.PredictMethod;
                    gproptshort.PredictMethod = PredictMethod;
                end

                if isfield(G, 'ActiveSetMethod')
                    ActiveSetMethod = G.ActiveSetMethod;
                    gproptshort.ActiveSetMethod = ActiveSetMethod;
                end

                if isfield(G, 'FitMethod')
                    FitMethod = G.FitMethod;
                    gproptshort.FitMethod = FitMethod;
                end

            end

            %Gaussian process regression
            if ~isempty(gpropts)
                gprMdl = fitrgp(X, y, gpropts{:});
            else
                gprMdl = fitrgp(X, y);
            end

            %compact the model
            cgprMdl = compact(gprMdl);

            %predictions ("interpolated" points)
            switch PredictMethod
                case 'fic'
                    [ypred, ysd, yint] = predict(gprMdl, X2);
                case 'bcd'
                    ypred = predict(cgprMdl, X2);
                otherwise
                    [ypred, ysd, yint] = predict(gprMdl, X2);
            end

            switch PredictMethod
                case 'fic'
                    mdlcmd = @(gprMdl, X2) predict(gprMdl, X2);
                    interpfn = @(qm2, nA2) interp_gpr(gprMdl, qm2, nA2, oref, projtol, usv, zeroQ);
                    mdlspec = var_names(gprMdl, cgprMdl, gpropts, ysd, yint);
                otherwise
                    mdlcmd = @(cgprMdl, X2) predict(cgprMdl, X2);
                    interpfn = @(qm2, nA2) interp_gpr(cgprMdl, qm2, nA2, oref, projtol, usv, zeroQ);
                    mdlspec = var_names(cgprMdl, gpropts, ysd, yint);
            end

            %model-specific parameters
            if exist('gproptshort', 'var') == 1

                if ~isempty(fieldnames(gproptshort))

                    if isfield(G, 'HyperparameterOptimizationOptions')
                        mdlparsspec = var_names(maxhyperobj, gproptshort);
                    else
                        mdlparsspec = var_names(gproptshort);
                    end

                else
                    evalc([(gproptstruct{end}) ' = G.(gproptstruct{end}']);
                    warning(['no tracked options were contained in user-specified gpropts. Consider tracking all. Tracking added via evalc() for ' ...
                            gproptstruct{end - 1}])

                    if isstruct(gproptshort.(gproptstruct{end}))
                        error(['Tracking automatically added for ' gproptstruct{end - 1} ' but the added tracking option cannot be a struct'])
                    end

                    mdlparsspec = struct(gproptstruct{end - 1}, gproptshort.gproptstruct{end});
                end

            else
                mdlparsspec = var_names(PredictMethod);
            end

            mdlspec.KernelInformation = cgprMdl.KernelInformation;
            mdlparsspec.KernelType = cgprMdl.KernelInformation.Name;
            mdlparsspec.KernelParameters = cgprMdl.KernelInformation.KernelParameters;
            mdlparsspec.KernelParameterNames = cgprMdl.KernelInformation.KernelParameterNames;
            mdlparsspec.Beta = cgprMdl.Beta;
            mdlparsspec.Sigma = cgprMdl.Sigma;

        case 'idw' % inverse distance weighting
            %whether to remove degenerate dimension or not
            if projQ
                X = ppts;
                X2 = ppts2;
            else
                X = o;
                X2 = o2;
            end

            r = nv.r;

            L = 2; %norm-power (i.e. L == 2 --> Euclidean norm)
            %different from Tovar's FEX idw.m implementation, but should be
            %similar or same output
            [ypred, W, r, nints, numnonints, int_fraction] = idw(X, X2, y, r, L);

            mdlcmd = @(X, X2, propList, r, L) idw(X, X2, propList, r, L);
            interpfn = @(qm2, nA2) interp_idw(X, qm2, nA2, y, r, L);

            %model-specific variables
            mdlspec = var_names(L, W, r, nints, numnonints, int_fraction);
            %model-specific parameters
            mdlparsspec = var_names(L, r, nints, numnonints, int_fraction);

        case 'nn' %nearest neighbors
            nnList = dsearchn(ppts, ppts2);

            mdlcmd = @(ppts, ppts2, propList) propList(dsearchn(ppts, ppts2));
            interpfn = @(qm2, nA2) interp_nn(ppts, qm2, nA2, projtol, usv, y);

            %assign NN property values
            ypred = y(nnList);

            %model-specific variables
            mdlspec.nnList = nnList;
            %model-specific parameters
            mdlparsspec = struct();

        case 'avg'
            % "interpolation" (just constant model)
            [ypred, yavg] = interp_avg(y, npredpts);

            mdlcmd = @(propList, ndatapts) interp_avg(propList, ndatapts);
            interpfn = @(qm2, nA2) repelem(yavg, npredpts, 1); %any new point gets assigned yavg

            %model-specific variables
            mdlspec.yavg = yavg;
            %model-specific parameters
            mdlparsspec = struct();
    end

    runtime = toc; %time elapsed to do the interpolation (method-specific portion)

    %% append extra general variables
    %parity variables
    parity.ypred = ypred;
    parity.ytrue = ytrue;

    %model
    mdlgen.ypred = ypred;
    mdlgen.mdlcmd = mdlcmd;
    mdlgen.interpfn = interpfn;
    mdlgen.runtime = runtime;
    mesh.ids = ids;
    mdlgen.mesh = mesh;
    mdlgen.data = data;
    mdlgen.parity = parity;
    %parameters
    mdlparsgen.runtime = runtime;
    mdlparsgen.parity = parity;

    %% concatenate structs
    %model variables
    mdl = structhorzcat(mdlgen, mdlspec);
    %model parameters
    mdlpars = structhorzcat(mdlparsgen, mdlparsspec);

end

function props = interp_gpr(gprMdl, qm2, nA2, oref, projtol, usv, zeroQ)
    % INTERP_GPR  interpolate using Gaussian Process Regression model (gprMdl) and predict()
    pts2 = get_pts(qm2, nA2, oref);
    ppts2 = get_ppts(pts2, projtol, usv, zeroQ);
    props = predict(gprMdl, ppts2);

end

function ppts = get_ppts(pts, projtol, usv, zeroQ)
    % GET_PPTS  get projected points via proj_down.m
    ppts = proj_down(pts, projtol, usv, 'zero', zeroQ);
end

function pts = get_pts(qm, nA, oref)
    % GET_PTS  get symmetrized octonions from qm/nA pairs
    pts = normr(get_octpairs(five2oct(qm, nA), 'oref', oref));
end

function [dataprops, facetprops, NNextrapID, nnList] = ...
        interp_bary_fast(ppts, ppts2, meshprops, databary, facetIDs)

    arguments
        ppts double
        ppts2 double
        meshprops double
        databary double
        facetIDs double
    end

    % INTERP_BARY_FAST  short-circuit barycentric interpolation (same input/prediction points, new mesh property values)
    % Cannot be used with new prediction points.

    %find NaN values & replace with NN values (NN extrapolation)
    [NNextrapID, ~] = isnan(databary);
    nnList = dsearchn(ppts2(NNextrapID), ppts);
    d = size(databary, 2);

    %properties of each vertex of each intersecting facet
    facetprops = get_facetprops(ppts, meshprops, facetIDs);

    % adjust NaN values to do NN extrapolation
    % e.g. databary == [NaN NaN NaN], facetprops == [NaN NaN NaN]
    % --> [1 0 0], [1.213 0 0], such that dot([1 0 0],[1.213 0 0])
    % == 1.213, where 1.213 is the extrapolated NN value
    databary(NNextrapID, 1) = 1;
    databary(NNextrapID, 2:d) = 0;
    facetprops(NNextrapID, 1) = meshprops(nnList);
    facetprops(NNextrapID, 2:d) = 0;

    dataprops = dot(databary, facetprops, 2);

end

%% HELPER FUNCTION(S)
function facetprops = get_facetprops(ppts, props, facetIDs)

    arguments
        ppts double
        props double
        facetIDs double
    end

    ndatapts = size(facetIDs, 1);
    facetprops = nan([ndatapts size(ppts, 2)]);

    for i = 1:ndatapts
        vtxIDs = facetIDs(i, :);
        facetprops(i, :) = props(vtxIDs).'; %properties of vertices of facet
    end

end

function props = interp_bary(mesh, intfacetIDs, qm2, nA2, usv, zeroQ, barytype, barytol, projtol, nnMax, brkQ, NV)

    arguments
        mesh(1, 1) struct {mustContainFields(mesh, {'pts', 'ppts', 'props', 'sphK'})}
        intfacetIDs cell
        qm2(:, 4) double
        nA2(:, 3) double
        usv(1, 1) struct
        zeroQ(1, 1) logical = false
        barytype char {mustBeMember(barytype, {'spherical', 'planar'})} = 'spherical'
        barytol(1, 1) double = 0.2
        projtol(1, 1) double = 1e-4
        nnMax(1, 1) double = 10
        brkQ(1, 1) logical = false
        NV.databary = []
        NV.propList = []
    end

    % INTERP_BARY  interpolate using spherical or planar barycentric coordinates
    pts2 = get_pts(qm2, nA2);
    ppts2 = get_ppts(pts2, projtol, usv, zeroQ);
    %% when barycentric coordinates aren't specified (new predict points)
    if isempty(intfacetIDs)
        intfacetIDs = intersect_facet(mesh.ppts, mesh.sphK, ppts2, inttol, 'inttype', 'planar', 'nnMax', nnMax);
    end

    data = struct('pts', pts2, 'ppts', ppts2);

    if brkQ
        five = struct('q', qm2, 'nA', nA2);
        data.props = GB5DOF_setup(five);
    else
        data.props = nan(size(pts2));
    end

    props = get_interp(mesh, data, intfacetIDs, barytype, barytol);

end

function yq = interp_idw(X, qm2, nA2, y, r, L)

    arguments
        X(:, 8) double %input points
        qm2(:, 4) double %query misorientations
        nA2(:, 3) double %query BP normals
        y(:, 1) double %property values
        r double = [] %radius
        L double = 2 % default is Euclidean norm
    end

    % INTERP_IDW interpolate using inverse-distance weighting and misorientation/BP normal query input pairs
    Xq = get_pts(qm2, nA2);
    yq = idw(X, Xq, y, r, L);

end

function props = interp_nn(ppts, qm2, nA2, projtol, usv, zeroQ, propList)
    % INTERP_NN  nearest neighbor interpolation
    if isempty(usv)
        error('usv should not be empty')
    end

    pts = get_pts(qm2, nA2);
    props = propList(dsearchn(ppts, get_ppts(pts, projtol, usv, zeroQ)));

end

function [gitcommit, comment, warnedQ] = get_gitcommit()
    % GET_GITCOMMIT  get git commit version (or return empty '' if error)
    [status, cmdout] = system('git rev-parse HEAD');

    if status == 0
        gitcommit = cmdout(1:7);
        comment = cmdout(9:end);
    else
        gitcommit = '';
        comment = '';
    end

    [status, cmdout] = system('git status --untracked-files=no --porcelain');

    if status == 0

        if ~isempty(cmdout)
            warning('Working directory not clean (i.e. uncommitted/unpushed) files exist. Use !git commit -am "<message>", then !git push')
            warnedQ = true;
        else
            warnedQ = false;
        end

    end

end

function uuid = get_uuid()
    % GET_UUID  get unique ID (8 characters, mixture of numbers and letters) via java.util.UUID.randomUUID
    temp = java.util.UUID.randomUUID;
    uuidtmp = char(temp.toString);
    uuid = uuidtmp(1:8);
end

function Sout = structhorzcat(S)

    arguments (Repeating)
        S(1, 1) struct
    end

    % STRUCTHORZCAT  "horizontally" concatenate structures with different variables
    %  (no overlapping variables, scalar struct).
    %--------------------------------------------------------------------------
    % Author(s): Sterling Baird
    %
    % Date: 2020-09-05
    %
    % Inputs:
    %  S - struct, where each struct corresponds to a set of e.g.
    %  conditions/parameters/results for a unique experiment, and each struct
    %  can have different numbers of rows and same and/or different variables*
    %
    % Outputs:
    %  Sout - "horizontally" concatenated structure
    %
    % Usage:
    %  Sout = tblvertcat(S1,S2);
    %  Sout = tblvertcat(S1,S2,S3);
    %
    % Dependencies:
    %
    % Notes:
    %--------------------------------------------------------------------------

    % convert structures to tables
    T = cellfun(@(S) struct2table(S, 'AsArray', true), S, 'UniformOutput', false);
    % concatenate tables
    Tcat = horzcat(T{:});
    % convert back to structure
    Sout = table2struct(Tcat);
end

function out = var_names(varargin)
    % VAR_NAMES  take variables and output a combined struct with each of the variable names
    %--------------------------------------------------------------------------
    % https://www.mathworks.com/matlabcentral/answers/79281#answer_89015
    %--------------------------------------------------------------------------
    for n = 1:nargin
        out.(inputname(n)) = varargin{n};
    end

end

function rad = deg2rad(deg)
    rad = deg / 180 * pi;
end

function A = allcomb(varargin)
    % ALLCOMB  All combinations (Cartesian Product)
    %    B = ALLCOMB(A1,A2,A3,...,AN) returns all combinations of the elements
    %    in the arrays A1, A2, ..., and AN. B is P-by-N matrix where P is the product
    %    of the number of elements of the N inputs.
    %    This functionality is also known as the Cartesian Product. The
    %    arguments can be numerical and/or characters, or they can be cell arrays.
    %
    %    Examples:
    %       allcomb([1 3 5],[-3 8],[0 1]) % numerical input:
    %       % -> [ 1  -3   0
    %       %      1  -3   1
    %       %      1   8   0
    %       %        ...
    %       %      5  -3   1
    %       %      5   8   1 ] ; % a 12-by-3 array
    %
    %       allcomb('abc','XY') % character arrays
    %       % -> [ aX ; aY ; bX ; bY ; cX ; cY] % a 6-by-2 character array
    %
    %       allcomb('xy',[65 66]) % a combination -> character output
    %       % -> ['xA' ; 'xB' ; 'yA' ; 'yB'] % a 4-by-2 character array
    %
    %       allcomb({'hello','Bye'},{'Joe', 10:12},{99999 []}) % all cell arrays
    %       % -> {  'hello'  'Joe'        [99999]
    %       %       'hello'  'Joe'             []
    %       %       'hello'  [1x3 double] [99999]
    %       %       'hello'  [1x3 double]      []
    %       %       'Bye'    'Joe'        [99999]
    %       %       'Bye'    'Joe'             []
    %       %       'Bye'    [1x3 double] [99999]
    %       %       'Bye'    [1x3 double]      [] } ; % a 8-by-3 cell array
    %
    %    ALLCOMB(..., 'matlab') causes the first column to change fastest which
    %    is consistent with matlab indexing. Example:
    %      allcomb(1:2,3:4,5:6,'matlab')
    %      % -> [ 1 3 5 ; 1 4 5 ; 1 3 6 ; ... ; 2 4 6 ]
    %
    %    If one of the N arguments is empty, ALLCOMB returns a 0-by-N empty array.
    %
    %    See also NCHOOSEK, PERMS, NDGRID
    %         and NCHOOSE, COMBN, KTHCOMBN (Matlab Central FEX)

    % Tested in Matlab R2015a and up
    % version 4.2 (apr 2018)
    % (c) Jos van der Geest
    % email: samelinoa@gmail.com

    % History
    % 1.1 (feb 2006), removed minor bug when entering empty cell arrays;
    %     added option to let the first input run fastest (suggestion by JD)
    % 1.2 (jan 2010), using ii as an index on the left-hand for the multiple
    %     output by NDGRID. Thanks to Jan Simon, for showing this little trick
    % 2.0 (dec 2010). Bruno Luong convinced me that an empty input should
    % return an empty output.
    % 2.1 (feb 2011). A cell as input argument caused the check on the last
    %      argument (specifying the order) to crash.
    % 2.2 (jan 2012). removed a superfluous line of code (ischar(..))
    % 3.0 (may 2012) removed check for doubles so character arrays are accepted
    % 4.0 (feb 2014) added support for cell arrays
    % 4.1 (feb 2016) fixed error for cell array input with last argument being
    %     'matlab'. Thanks to Richard for pointing this out.
    % 4.2 (apr 2018) fixed some grammar mistakes in the help and comments

    narginchk(1, Inf);
    NC = nargin;

    % check if we should flip the order
    if ischar(varargin{end}) && (strcmpi(varargin{end}, 'matlab') || strcmpi(varargin{end}, 'john'))
        % based on a suggestion by JD on the FEX
        NC = NC - 1;
        ii = 1:NC; % now first argument will change fastest
    else
        % default: enter arguments backwards, so last one (AN) is changing fastest
        ii = NC:-1:1;
    end

    args = varargin(1:NC);

    if any(cellfun('isempty', args)) % check for empty inputs
        warning('ALLCOMB:EmptyInput', 'One of more empty inputs result in an empty output.');
        A = zeros(0, NC);
    elseif NC == 0 % no inputs
        A = zeros(0, 0);
    elseif NC == 1 % a single input, nothing to combine
        A = args{1}(:);
    else
        isCellInput = cellfun(@iscell, args);

        if any(isCellInput)

            if ~all(isCellInput)
                error('ALLCOMB:InvalidCellInput', ...
                'For cell input, all arguments should be cell arrays.');
            end

            % for cell input, we use to indices to get all combinations
            ix = cellfun(@(c) 1:numel(c), args, 'un', 0);

            % flip using ii if last column is changing fastest
            [ix{ii}] = ndgrid(ix{ii});

            A = cell(numel(ix{1}), NC); % pre-allocate the output

            for k = 1:NC
                % combine
                A(:, k) = reshape(args{k}(ix{k}), [], 1);
            end

        else
            % non-cell input, assuming all numerical values or strings
            % flip using ii if last column is changing fastest
            [A{ii}] = ndgrid(args{ii});
            % concatenate
            A = reshape(cat(NC + 1, A{:}), [], NC);
        end

    end

    %% COPYRIGHT
    % Copyright (c) 2018, Jos (10584)
    % All rights reserved.
    %
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are
    % met:
    %
    %     * Redistributions of source code must retain the above copyright
    %       notice, this list of conditions and the following disclaimer.
    %     * Redistributions in binary form must reproduce the above copyright
    %       notice, this list of conditions and the following disclaimer in
    %       the documentation and/or other materials provided with the distribution
    %
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    % POSSIBILITY OF SUCH DAMAGE.

end

function o = five2oct(qm, nA, epsijk)

    arguments
        qm(:, 4) double {mustBeFinite, mustBeReal}
        nA(:, 3) double {mustBeFinite, mustBeReal}
        epsijk(1, 1) double = 1
    end

    % FIVE2OCT Convert 5DOF coordinates (qm,nA) to octonions (o).
    % Misorientation quaternion and grain A crystal frame boundary plane normal
    % pointing outward from grain A towards grain B to octonion with BP normal
    % = [0 0 1];
    %--------------------------------------------------------------------------
    % Author(s): Sterling Baird
    % Date: 2020-12-05
    %
    % Inputs:
    %  qm - Misorientation quaternion
    %    resp.
    %  nA - boundary plane normal (grain A crystal frame) pointing from grain
    %  A to grain B
    %
    % Outputs:
    %   o - octonion, with BP normal = [0 0 1]
    %
    % Usage:
    %  o = five2oct(qm,nA)
    %
    % Dependencies:
    %  qinv.m
    %  qmA2oct.m
    %
    % References:
    %  (1) Francis, T.; Chesser, I.; Singh, S.; Holm, E. A.; De Graef, M. A
    %  Geodesic Octonion Metric for Grain Boundaries. Acta Materialia 2019,
    %  166, 135–147. https://doi.org/10.1016/j.actamat.2018.12.034.

    %--------------------------------------------------------------------------

    npts = size(qm, 1);
    % assign identity quaternion to grain A
    pA = repmat([1 0 0 0], npts, 1);
    % assign misorientation quaternion to grain B
    pB = qm;
    % because pA is identity quaternion, mA == nA
    mA = nA;

    if epsijk == -1
        %correct for misorientation convention
        pA = qinv(pA);
        pB = qinv(pB);
    end

    % conversion
    o = qmA2oct(pA, pB, nA, epsijk);

end

function propList = GB5DOF_setup(pA, pB, mA, mat, epsijk, nv)

    arguments
        pA(:, 4) = []
        pB(:, 4) = []
        mA(:, 3) = [0 0 1] %default is octonion convention
        mat char = 'Ni'
        epsijk(1, 1) double = 1
        nv.o double = []
    end

    %GB5DOF_SETUP  Compute 5DOF GB energy from BRK function
    %--------------------------------------------------------------------------
    % Author(s): Oliver Johnson, Sterling Baird
    % Date: 2020-07-27
    %
    % Inputs:
    %  five - struct containing at minimum misorientation quaternions (q) and
    %  boundary plane normals in crystal reference frame of grain A (nA)
    %  pointing towards grain B.
    %
    % Outputs:
    %  propList - grain boundary energies computed at each grain boundary in
    %  five
    %
    % Usage:
    %  propList = GB5DOF_setup(five)
    %
    % Dependencies:
    %  constructGBMatrices.m
    %  GB5DOF.m
    %
    % Notes:
    %  Find a good verification/ground truth (i.e. here is a GBE for a GB with
    %  this representation.)
    %--------------------------------------------------------------------------

    % Compute GB matrices
    % if ~isempty(five)
    %     pB = vertcat(five.q);
    %     pA = repmat([1 0 0 0],nGB,1);
    %     mA = vertcat(five.nA).';
    % end

    o = nv.o;

    if isempty([pA; pB])
        assert(~isempty(o), 'specify either pA and pB or o')
        pA = o(:, 1:4);
        pB = o(:, 5:8);
    end

    pA = normr(pA);
    pB = normr(pB);

    npts = size(pB, 1);

    if isempty(mA)
        mA = [0 0 1];
    end

    if size(mA, 1) == 1
        mA = repmat(mA, npts, 1);
    end

    if isempty(pA)
        % case that pB is actually misorientation
        pA = repmat([1 0 0 0], npts, 1);
    end

    if epsijk == -1
        pA = qinv(pA);
        pB = qinv(pB);
    end

    [omA, omB] = deal(zeros(3, 3, npts));

    parfor i = 1:npts
        mAtmp = mA(i, :);
        R = vecpair2rmat(mAtmp, [1 0 0], 1);
        qR = om2qu(R, 1);
        pAtmp = pA(i, :);
        pBtmp = pB(i, :);
        qA = qmult(qR, pAtmp, 1);
        qB = qmult(qR, pBtmp, 1);
        omA(:, :, i) = qu2om(qA, 1);
        omB(:, :, i) = qu2om(qB, 1);
    end

    %Calculate GB Energies
    % mat = 'Ni'; %'Al', 'Au', 'Cu'
    E(npts) = struct;
    E(1).(mat) = [];

    parfor k = 1:npts
        E(k).(mat) = GB5DOF(omA(:, :, k), omB(:, :, k), mat);
    end

    propList = vertcat(E.(mat));
end

function en = GB5DOF(P, Q, AlCuParameter, eRGB)
    %  GB5DOF  computes the energy of an arbitrary boundary in FCC metals (BRK energy function)
    %
    %     en = GB5DOF(P,Q,AlCuParameter) computes the energy of a boundary
    %     described by two rotation matrices P and Q. The character string
    %     variable AlCuParameter can take on one of four values: 'Al', 'Ni',
    %     'Au', or 'Cu'.  The output variable en is the computed boundary energy
    %     in units J/m^2.
    %
    %     en = GB5DOF(P,Q,AlCuParameter,eRGB) calculates the energy of a boundary
    %     in a hypothetical FCC metal defined by the numerical values of input
    %     parameters AlCuParameter and eRGB.  eRGB defines the scale of boundary
    %     energy variations in the hypothetical metal and should be given in
    %     units J/m^2. When eRGB is defined, parameter AlCuParameter should
    %     take on a numerical value ranging from 0.0 (for aluminum) to 1.0
    %     (for copper).
    %
    %     P and Q are two properly normalized 3x3 rotation matrices defining
    %     orientations of the two grains with respect to a laboratory (sample)
    %     frame. For any vector V expressed in the cube frame of grain P (or Q),
    %     P*V (or Q*V) expresses the same vector in the laboratory frame. By
    %     convention, the first row of P (or Q) is the row vector of the boundary
    %     plane normal Np = P(1,:) (or Nq = Q(1,:)) written in the cube frame
    %     of grain P (or Q). Thus, P*Np' = Q*Nq' = [1 0 0]'.
    %
    %     Examples
    %
    %       With P and Q matrices defined as follows
    %
    %           P = [ 0.5774    0.5774    0.5774 ;
    %                 0.7071   -0.7071         0 ;
    %                 0.4082    0.4082   -0.8165 ]
    %
    %           Q = [ 0.5774    0.5774    0.5774 ;
    %               [-0.7071    0.7071         0 ;
    %                -0.4082   -0.4082    0.8165 ]
    %
    %       en = GB5DOF(P,Q,'Ni') returns en = 0.0624 which is the energy in
    %       metal Ni of the coherent twin boundary defined by matrices P and Q.
    %
    %       With the same matrices P and Q, en = GB5DOF(P,Q,0.768,1.445) returns
    %       the same value en = 0.0624. In this example numerical parameters
    %       AlCuParameter = 0.768 and eRGB = 1.445 have values exactly matching
    %       the best fit values of the same parameters for FCC metal Ni.

    geom100 = distances_to_set(P, Q, '100'); % Generate geometry parameters
    geom110 = distances_to_set(P, Q, '110');
    geom111 = distances_to_set(P, Q, '111');

    if ~exist('eRGB', 'var') || isempty(eRGB)
        parvec = makeparvec(AlCuParameter); % Option 2
    else
        parvec = makeparvec(AlCuParameter, eRGB); % Option 1
    end

    en = weightedmeanenergy(geom100, geom110, geom111, parvec); % Calculate the energy
end

function geom = distances_to_set(P, Q, whichaxes, dismax)
    % geom = distances_to_set(P,Q,whichaxes,dismax)
    %
    % Calculates the geometry parameters for a given grain boundary relative to
    % a given set of axes.
    %
    % P and Q are rotation matrices giving the orientations of the two grains.
    % The grain boundary normal is fixed at [1,0,0].
    %
    % whichaxes is one of '100', '110', or '111'
    %
    % dismax is an optional parameter specifying the maximum distance to
    % include. It defaults to slightly less than 1, which is the largest
    % allowable distance before some anomalies start to appear.
    %
    % Result geom is a 4xn matrix, where the rows are distance, ksi, eta, and
    % phi. It keeps all of the hits up to a distance of dismax.
    %
    % distance is 2*sin(delta/2) where delta is the angle of closest approach
    % between a misorientation axis and one of the axes.  Note there are 24
    % representations of the rotations and 3, 6, or 4 equivalent high-symmetry
    % axes, so it calculates as many as 144 distances.  But only ones below
    % the cutoff dismax are kept.
    %
    % Once it's picked the closest approximation to the boundary for a given
    % axis and coset element, it finds the parameters ksi, eta, phi defining
    % that idealized boundary (since the axis is defined, it's a 3-space).
    %
    % These are:
    % phi, the angle between the rotation axis and the boundary plane normal
    % (taken as the mean of the normals represented in the two actual grain
    % orientations, which works when dismax is less than 1)
    %
    % ksi, the misorientation angle
    %
    % eta, a parameter giving the second axis of the boundary plane normal in
    % terms of specified directions ('dirs') perpendicular to each
    % high-symmetry axis.

    if ~exist('dismax', 'var') || isempty(dismax)
        dismax = 0.999999; % Force the distance to be strictly less than one, allowing for roundoff
    end % Note if dismax >= 1, you're likely to get a warning about m1 being singular.

    switch whichaxes
        case {'110'}
            % Define 110 axes, normalize
            axes = [1 1 1 1 0 0;
                1 -1 0 0 1 1;
                0 0 1 -1 1 -1] / sqrt(2);

            % Define a crystal direction perpendicular to each rotation axis.
            % The formalism demands that this be an axis of at least two-fold
            % symmetry.
            dirs = [0 0 0 0 1 1;
                0 0 1 1 0 0;
                1 1 0 0 0 0];
        case {'111'}

            % Define 111 axes, normalize
            axes = [1 1 -1 -1;
                1 -1 1 -1;
                1 -1 -1 1] / sqrt(3);

            dirs = [1 1 1 1;
                -1 1 1 -1;
                0 0 0 0] / sqrt(2);

        case {'100'}
            % Define 100 axes, normalize
            axes = [1 0 0;
                0 1 0;
                0 0 1];

            dirs = [0 0 1;
                1 0 0;
                0 1 0];
        otherwise
            error('Undefined axis set')
    end

    naxes = size(axes, 2);
    period = pi * naxes / 6;

    %  Define the symmetry operators

    rotX90 = [1 0 0; %  Rotation by +90 degrees around X axis
        0 0 -1;
        0 1 0];

    rotY90 = [0 0 1; %  Rotation by +90 degrees around Y axis
        0 1 0;
        -1 0 0];

    rotZ90 = [0 -1 0; %  Rotation by +90 degrees around Z axis
        1 0 0;
        0 0 1];

    rotZ90m = [0 1 0; %  Rotation by -90 degrees around Z axis
            -1 0 0;
            0 0 1];

    % Create 24 symmetry equivalent variants of Q
    % This is the coset appropriate for the rotation convention where Q'*P
    % is the misorientation represented in the grain frame.  If you're
    % getting odd results, e.g. misorientations that you know are CSL are
    % coming out entirely wrong, you may be using the opposite convention;
    % try replacing P and Q with P' and Q'.
    V = cell(24, 1);

    V{1} = Q;
    V{2} = V{1} * rotX90; % Rotate the vectors three times around X by +90 degrees
    V{3} = V{2} * rotX90;
    V{4} = V{3} * rotX90;

    for j = 1:12 % Rotate three times around Y by +90 degrees
        V{j + 4} = V{j} * rotY90;
    end

    for j = 1:4
        V{j + 16} = V{j} * rotZ90; % Rotate three times around Z by +90 degrees
        V{j + 20} = V{j} * rotZ90m; % Rotate three times around Z by -90 degrees
    end

    % Preallocate all parameter lists at their maximum possible sizes.
    % Redundant representations will be removed at the end.
    distances = zeros(1, 24 * naxes);
    phis = zeros(1, 24 * naxes);
    ksis = zeros(1, 24 * naxes);
    etas = zeros(1, 24 * naxes);

    thisindex = 0; % Number of hits found so far

    % Step through all combinations of symmetrically-equivalent axes and coset
    % elements V{j}.
    for i = 1:naxes

        ax = axes(:, i); % ax is a high-symmetry axis

        dir = dirs(:, i); %  This is the pivot vector used to partition
        %  the rotation around axis "i"
        dir2 = cross(ax, dir); %  Completing the orthonormal coordinate set.
        %  theta1 and theta2 are defined in the plane
        %  spanned by (dir,dir2).

        for j = 1:24 % For each symmetry-related variant of the second grain

            Q = V{j};
            R = Q' * P; %  This rotates any vector in cube P into a vector in cube Q

            % Edited by SGB, 2020-07-09, to prevent imaginary output later
            %R		= round(R,12);
            R(abs(R) < 1e-12) = 0;
            R(abs(R - 1) < 1e-12) = 1;
            R(abs(R + 1) < 1e-12) = -1;

            q = mat2quat(R); % Calculation from here on out is much easier with quaternions.
            axi = q(2:4)' / sqrt(sum(q(2:4).^2)); % Normalized rotation axis
            psi = 2 * acos(q(1)); % Rotation angle

            dotp = axi * ax;

            % Compute rotational distance from boundary P/Q to the rotation set "i"
            % This formula produces 2*sin(delta/2), where delta is the angle of
            % closest approach.
            dis = 2 * sqrt(abs(1 - dotp * dotp)) * sin(psi / 2);

            if dis < dismax
                thisindex = thisindex + 1;

                theta = 2 * atan(dotp * tan(psi / 2)); % angle of rotation about ax that most closely approximates R

                % Compute the normal of the best-fitting GB in grain 1
                n1 = P(1, :)';
                n2 = Q(1, :)';

                RA = quat2mat([cos(theta / 2); sin(theta / 2) * ax]);
                % RA is the rotation about ax that most closely approximates R

                % From this point on we're dealing with the idealized rotation RA, not
                % the original rotation R.
                m1 = n1 + RA' * n2;

                % The next problem comes up only for very large distances,
                % which are normally cut off
                if norm(m1) < 0.000001
                    disp('m1 is singular!!!')
                end

                m1 = m1 / norm(m1); % Halfway between the two normal vectors from the two grains
                m2 = RA * m1; % And the same represented in the other grain

                % Compute the inclination angle for the common rotation axis
                phi = real(acos(abs(m1' * ax))); % "real" because of numerical problems when they're exactly parallel

                % Partition the total rotation angle "theta"
                if abs(ax' * m1) > 0.9999 % Check if the best-fitting GB is pure twist
                    theta1 =- theta / 2; % eta is meaningless for a twist boundary.
                    theta2 = theta / 2;
                else

                    try
                        theta1 = atan2(dir2' * m1, dir' * m1);
                        theta2 = atan2(dir2' * m2, dir' * m2);
                    catch
                        disp('')
                    end

                    % It's projecting m1 and m2 into the plane normal to ax and
                    % then determining the rotation angles of them relative to
                    % dir.

                end

                % Reduce both angles to interval (-period/2,period/2],
                % semi-open with a small numerical error.
                theta2 = theta2 - round(theta2 / period) * period;
                theta1 = theta1 - round(theta1 / period) * period;

                % This implements the semi-open interval in order to avoid an
                % annoying numerical problem where certain representations are
                % double-counted.
                if abs(theta2 + period / 2) < 0.000001
                    theta2 = theta2 + period;
                end

                if abs(theta1 + period / 2) < 0.000001
                    theta1 = theta1 + period;
                end

                % Since this is only being run on fcc elements, which are
                % centrosymmetric, and all dir vectors are 2-fold axes, then
                % the operations of swapping theta1 and theta2, and of
                % multilying both by -1, are symmetries for the energy
                % function. This lets us fold everything into a small right
                % triangle in (ksi,eta) space:
                ksi = abs(theta2 - theta1);
                eta = abs(theta2 + theta1);

                % And store them in the vectors
                distances(thisindex) = dis;
                ksis(thisindex) = ksi;
                etas(thisindex) = eta;
                phis(thisindex) = phi;
            end

        end

    end

    % Dump the excess pre-allocated ones and sort the rest in order of distance
    [distances, sortindex] = sort(distances(1:thisindex));
    ksis = ksis(sortindex);
    etas = etas(sortindex);
    phis = phis(sortindex);

    % Clean up redundancy.  Double-counting the same representation of one
    % boundary messes up the weighting functions in weightedmeanenergy.m

    % First round everything to 1e-6, so that negligible numerical
    % differences are dropped
    distances = 1e-6 * round(distances * 1e6);
    ksis = 1e-6 * round(ksis * 1e6);
    etas = 1e-6 * round(etas * 1e6);
    phis = 1e-6 * round(phis * 1e6);

    % And finally create the 4 x thisindex array of geometrical parameters
    geom = unique([distances', ksis', etas', phis'], 'rows')';

end

function q = mat2quat(m)
    % q = mat2quat(m)
    %
    % Auxiliary function converts a rotation matrix, assumed orthonormal, into
    % a unit quaternion.
    t = m(1, 1) + m(2, 2) + m(3, 3);
    e0 = sqrt(1 + t) / 2;

    if t > -0.999999999
        e = [m(2, 3) - m(3, 2); m(3, 1) - m(1, 3); m(1, 2) - m(2, 1)] / (4 * e0);
    else
        e0 = 0;
        e3 = sqrt(- (m(1, 1) + m(2, 2)) / 2);

        if abs(e3) > 2e-8 % Check for singularity, allowing numerical error
            e = [m(1, 3) / (2 * e3); m(2, 3) / (2 * e3); e3];
        else
            e1 = sqrt((m(1, 1) + 1) / 2);

            if e1 ~= 0
                e = [e1; m(2, 1) / (2 * e1); 0];
            else
                e = [0; 1; 0];
            end

        end

    end

    q = [e0; -e];
end

function m = quat2mat(q)
    % m = quat2mat(q)
    %
    % Auxiliary function converts a quaternion into a rotation matrix with no
    % assumption about normalization.
    e0 = q(1);
    e1 = q(2);
    e2 = q(3);
    e3 = q(4);

    m = [e0^2 + e1^2 - e2^2 - e3^2, 2 * (e1 * e2 - e0 * e3), 2 * (e1 * e3 + e0 * e2); ...
            2 * (e1 * e2 + e0 * e3), e0^2 - e1^2 + e2^2 - e3^2, 2 * (e2 * e3 - e0 * e1); ...
            2 * (e1 * e3 - e0 * e2), 2 * (e2 * e3 + e0 * e1), e0^2 - e1^2 - e2^2 + e3^2] ...
        / (e0^2 + e1^2 + e2^2 + e3^2);
end

function [par43, AlCuparameter] = makeparvec(AlCuparameter, eRGB, par42Al, par42Cu)
    % [par43,AlCuparameter] = makeparvec(AlCuparameter,eRGB,par42Al,par42Cu)
    %
    % Creates a 43-parameter vector as used by weightedmeanenergy
    %
    % Arguments are:
    % AlCuparameter:  Position on the Al-Cu axis, where 0 is Al and 1 is Cu
    % (this parameter is capital Phi in the journal article).
    % This is related to eSF/eRGB, where eSF is the stacking-fault energy.
    % Optionally, AlCu parameter is a string 'Al', 'Ni', 'Au', or 'Cu', which
    % then sets all other parameters automatically.  You can call it with
    % just this parameter if you wish.
    %
    % eRGB:  Energy of a "random" grain boundary in J/m^2
    %
    % There are some additional options that are still written into the
    % function but are currently not used:
    % par42Al:  The 42 dimensionless parameters for Al
    %
    % par42Cu:  The 42 dimensionless parameters for Cu
    %
    % Note a majority of the entries for par42Al and par42Cu are normally
    % identical.
    %
    % All parameters have default values.  Defaults for par42Al and par42Cu are
    % the values found by numerical fitting to the 388*4 boundaries.
    % eRGB and AlCuparameter default to the values for Cu.
    %
    % Optionally returns the numerical AlCuparameter so the user can read the
    % default value for each element.

    if ~exist('eRGB', 'var') || isempty(eRGB)
        eRGB = 1.03669431227427; % Value for Cu
    end

    if ~exist('AlCuparameter', 'var') || isempty(AlCuparameter)
        AlCuparameter = 1; % Value for Cu
    end

    if ~exist('par42Al', 'var') || isempty(par42Al)
        par42Al = [0.405204179289160; 0.738862004021890; 0.351631012630026; 2.40065811939667; 1.34694439281655; 0.352260396651516; 0.602137375062785; 1.58082498976078; 0.596442399566661; 1.30981422643602; 3.21443408257354; 0.893016409093743; 0.835332505166333; 0.933176738717594; 0.896076948651935; 0.775053293192055; 0.391719619979054; 0.782601780600192; 0.678572601273508; 1.14716256515278; 0.529386201144101; 0.909044736601838; 0.664018011430602; 0.597206897283586; 0.200371750006251; 0.826325891814124; 0.111228512469435; 0.664039563157148; 0.241537262980083; 0.736315075146365; 0.514591177241156; 1.73804335876546; 3.04687038671309; 1.48989831680317; 0.664965104218438; 0.495035051289975; 0.495402996460658; 0.468878130180681; 0.836548944799803; 0.619285521065571; 0.844685390948170; 1.02295427618256];
    end

    if ~exist('par42Cu', 'var') || isempty(par42Cu)
        par42Cu = [0.405204179289160; 0.738862004021890; 0.351631012630026; 2.40065811939667; 1.34694439281655; 3.37892632736175; 0.602137375062785; 1.58082498976078; 0.710489498577995; 0.737834049784765; 3.21443408257354; 0.893016409093743; 0.835332505166333; 0.933176738717594; 0.896076948651935; 0.775053293192055; 0.509781056492307; 0.782601780600192; 0.762160812499734; 1.10473084066580; 0.529386201144101; 0.909044736601838; 0.664018011430602; 0.597206897283586; 0.200371750006251; 0.826325891814124; 0.0226010533470218; 0.664039563157148; 0.297920289861751; 0.666383447163744; 0.514591177241156; 1.73804335876546; 2.69805148576400; 1.95956771207484; 0.948894352912787; 0.495035051289975; 0.301975031994664; 0.574050577702240; 0.836548944799803; 0.619285521065571; 0.844685390948170; 0.0491040633104212];
    end

    if ischar(AlCuparameter)

        switch AlCuparameter
            case 'Ni'
                eRGB = 1.44532834613925;
                AlCuparameter = 0.767911805073948;
            case 'Al'
                eRGB = 0.547128733614891;
                AlCuparameter = 0;
            case 'Au'
                eRGB = 0.529912885175204;
                AlCuparameter = 0.784289766313152;
            case 'Cu'
                eRGB = 1.03669431227427;
                AlCuparameter = 1;
            otherwise
                error('Undefined element')
        end

    end

    par43 = [eRGB; (par42Al + AlCuparameter * (par42Cu - par42Al))];
end

function en = weightedmeanenergy(geom100, geom110, geom111, pars)
    % en = weightedmeanenergy(geom100,geom110,geom111,pars)
    % Calculate the energy for a single grain boundary.
    %
    % Input variables geom100, geom110, and geom111 are each matrices with 4
    % rows, giving the non-redundant representations of the boundary about each
    % set of axes as generated by distances_to_set.m.  See the comments in that
    % function for further information.  The rows are distance;ksi;eta;phi.
    %
    % Input variable pars is a 43-element vector as created by makeparvec.m.
    % This specifies all of the parameters needed for the 5DOF energy function
    % on a specific fcc metal.
    %
    % Return variable en is the energy in J/m^2.
    %
    % The physical relevance of the parameters is commented wherever they
    % appear, in this function and in the setxxx functions.

    % Pull out the parameters relevant to the weighting of the 100, 110, and
    % 111 sets
    eRGB = pars(1); % The only dimensioned parameter.  The energy scale, set by the energy of a random boundary.
    d0100 = pars(2); % Maximum distance for the 100 set.  Also the distance scale for the rsw weighting function.
    d0110 = pars(3); % Same for the 110 set
    d0111 = pars(4); % Same for the 111 set
    weight100 = pars(5); % Weight for the 100 set, relative to the random boundary
    weight110 = pars(6); % Same for 110
    weight111 = pars(7); % Same for 111

    offset = 0.00001; % Cutoff of weighting function at small d, purely for numerical purposes

    % The following three energy lists are in units of eRGB.
    e100 = set100(geom100, pars);
    e110 = set110(geom110, pars);
    e111 = set111(geom111, pars);

    d100 = geom100(1, :);
    d110 = geom110(1, :);
    d111 = geom111(1, :);

    % Now calculate the weights, in a manner designed to give an rsw-like
    % function of d.  Note it calculates a weight for every representation of
    % the boundary within each set.
    s100 = sin(pi / 2 * d100 / d0100);
    s100(d100 > d0100) = 1; % Weight saturates at zero above d0
    s100(d100 < offset * d0100) = offset * pi / 2; % Avoid calculation of NaN's, replace with something small but finite
    w100 = (1 ./ (s100 .* (1 - 0.5 * log(s100))) - 1) * weight100;

    s110 = sin(pi / 2 * d110 / d0110);
    s110(d110 > d0110) = 1;
    s110(d110 < offset * d0110) = offset * pi / 2;
    w110 = (1 ./ (s110 .* (1 - 0.5 * log(s110))) - 1) * weight110;

    s111 = sin(pi / 2 * d111 / d0111);
    s111(d111 > d0111) = 1;
    s111(d111 < offset * d0111) = offset * pi / 2;
    w111 = (1 ./ (s111 .* (1 - 0.5 * log(s111))) - 1) * weight111;

    en = eRGB * (sum(e100 .* w100) + sum(e110 .* w110) + sum(e111 .* w111) + 1) / (sum(w100) + sum(w110) + sum(w111) + 1);
end

function en = set100(geom100, pars)
    % en = set100(geom100,pars)
    %
    % Calculate the dimensionless contribution to the boundary based on the
    % nearby <100> rotations.  Meant to be called by weightedmeanenergy.m, but
    % also can be a stand-alone function for purposes of plotting cross
    % sections through the function.
    % Input variables geom100 and pars are as generated by distances_to_set.m
    % and makeparvec.m.  See comments in those functions for more information.

    pwr1 = pars(8); % 100 tilt/twist mix power law:  Twist
    pwr2 = pars(9); % 100 tilt/twist mix power law:  Tilt

    ksi = geom100(2, :);
    eta = geom100(3, :);
    phi = geom100(4, :);

    entwist = twist100(ksi, pars);
    entilt = atgb100(eta, ksi, pars);

    x = phi / (pi / 2);
    en = entwist .* (1 - x).^pwr1 + entilt .* x.^pwr2;

end

function en = twist100(ksi, pars)
    % en = twist100(ksi,pars)
    %
    % Dimensionless 100 twist contribution to the energy

    a = pars(10); % 100 twist maximum energy
    b = pars(10) * pars(11); % 100 twist rsw shape factor. The unusual split into two parameters is a holdover from an older version.

    perio = pi / 2; % the twist period
    ksi = mod(abs(ksi), perio); % rotation symmetry

    ksi(ksi > perio / 2) = perio - ksi(ksi > perio / 2);

    % Implement an rsw function of ksi
    sins = sin(2 * ksi);
    xlogx = sins .* log(sins);
    xlogx(isnan(xlogx)) = 0; % Force the limit to zero as x -> 0.
    en = a * sins - b * xlogx;

end

function en = atgb100(eta, ksi, pars)
    % en = atgb100(eta,ksi,pars)
    %
    % This function is a fit to the energies of all 100-tilt boundaries
    %

    pwr = pars(12); % 100 atgb interpolation power law

    period = pi / 2;

    en1 = stgb100(ksi, pars); % Value at eta = 0
    en2 = stgb100(period - ksi, pars); % Value at eta = pi/2

    % eta dependence is a power law that goes from the higher to the lower,
    % whichever direction that is
    select = en1 >= en2;
    en = zeros(size(ksi));

    en(select) = en1(select) - (en1(select) - en2(select)) .* (eta(select) / period).^pwr;
    en(~select) = en2(~select) - (en2(~select) - en1(~select)) .* (1 - eta(~select) / period).^pwr;
end

function en = stgb100(ksi, pars)
    % en = stgb100(ksi,pars)
    %
    % dimensionless 100 symmetric tilt energy
    %
    % This is implemented as a piecewise-rsw function, specified by energy
    % parameters en, angle breaks th, and shape factors a.

    en2 = pars(13); % peak before first Sigma5
    en3 = pars(14); % first Sigma5
    en4 = pars(15); % peak between Sigma5's
    en5 = pars(16); % second Sigma5
    en6 = pars(17); % Sigma17

    th2 = pars(18); % position of peak before first Sigma5
    th4 = pars(19); % position of peak between Sigma5's

    th6 = 2 * acos(5 / sqrt(34)); %Sigma17 rotation angle
    a12 = 0.5; % rsw shape factor.  In previous versions, these were allowed
    a23 = a12; % to vary, however there were too few vicinal boundaries in the
    a34 = a12; % ensemble to constrain them.  We found that forcing the great
    a45 = a12; % majority of them to be 0.5 helped to constrain the fit and
    a56 = a12; % produced reasonable results.  This holds true for most of the
    a67 = a12; % rsw shape factors throughout this code.
    %

    en1 = 0; % Sigma1 at left end
    en7 = 0; % Sigma1 at right end

    th1 = 0; % Sigma1 at left end
    th3 = acos(4/5); % first Sigma5
    th5 = acos(3/5); % second Sigma5
    th7 = pi / 2; % Sigma1 at right end
    %
    % And the rest is just the piecewise rsw function itself.
    en = zeros(size(ksi));
    select = (ksi <= th2);
    en(select) = en1 + (en2 - en1) * rsw(ksi(select), th1, th2, a12);
    select = (ksi >= th2 & ksi <= th3);
    en(select) = en3 + (en2 - en3) * rsw(ksi(select), th3, th2, a23);
    select = (ksi >= th3 & ksi <= th4);
    en(select) = en3 + (en4 - en3) * rsw(ksi(select), th3, th4, a34);
    select = (ksi >= th4 & ksi <= th5);
    en(select) = en5 + (en4 - en5) * rsw(ksi(select), th5, th4, a45);
    select = (ksi >= th5 & ksi <= th6);
    en(select) = en6 + (en5 - en6) * rsw(ksi(select), th6, th5, a56);
    select = (ksi >= th6 & ksi <= th7);
    en(select) = en7 + (en6 - en7) * rsw(ksi(select), th7, th6, a67);
end

function en = set110(geom110, pars)
    % en = set110(geom110,pars)
    %
    % Dimensionless contribution to energy from <110> rotations
    % Very similar to set100; see comments therein for general information.
    % Comments in this file will be limited to 110-specific information.

    pwr1 = pars(20); % 110 tilt/twist mix power law:  Twist
    pwr2 = pars(21); % 110 tilt/twist mix power law:  Tilt

    ksi = geom110(2, :);
    eta = geom110(3, :);
    phi = geom110(4, :);

    %
    entwist = twists110(ksi, pars);
    entilt = atgbs110(eta, ksi, pars);
    x = phi / (pi / 2);
    en = entwist .* (1 - x).^pwr1 + entilt .* x.^pwr2;
end

function en = atgbs110(eta, ksi, pars)
    % en = atgbs110(eta,ksi,pars)
    % See comments on set110.

    a = pars(26); % 110 atgb interpolation rsw shape factor
    %
    period = pi;
    en1 = stgbs110(ksi, pars);
    en2 = stgbs110(period - ksi, pars);

    en = zeros(size(eta));

    % Power-law interpolation did not work well in this case.  Did an rsw
    % function instead.
    select = en1 >= en2;

    en(select) = en2(select) + (en1(select) - en2(select)) .* rsw(eta(select), pi, 0, a);
    en(~select) = en1(~select) + (en2(~select) - en1(~select)) .* rsw(eta(~select), 0, pi, a);

end

function en = stgbs110(th, pars)
    % en = stgbs110(th,pars)
    % See comments on set110.

    en2 = pars(27); % peak between Sigma1 and Sigma3
    en3 = pars(28); % Coherent Sigma3 twin relative energy; one of the more important element-dependent parameters
    en4 = pars(29); % energy peak between Sigma3 and Sigma11
    en5 = pars(30); % Sigma11 energy
    en6 = pars(31); % energy peak between Sigma11 and Sigma1

    th2 = pars(32); % peak between Sigma1 and Sigma3
    th4 = pars(33); % peak between Sigma3 and Sigma11
    th6 = pars(34); % peak between Sigma11 and higher Sigma1

    a12 = 0.5;
    a23 = 0.5;
    a34 = 0.5;
    a45 = 0.5;
    a56 = 0.5;
    a67 = 0.5;
    %
    %
    en1 = 0;
    en7 = 0;

    th1 = 0;
    th3 = acos(1/3); % Sigma3
    th5 = acos(-7/11); % Sigma11
    th7 = pi;
    %
    th = pi - th; % This is a legacy of an earlier (ksi,eta) mapping
    %
    en = zeros(size(th));

    select = th <= th2;
    en(select) = en1 + (en2 - en1) .* rsw(th(select), th1, th2, a12);

    select = th >= th2 & th <= th3;
    en(select) = en3 + (en2 - en3) .* rsw(th(select), th3, th2, a23);

    select = th >= th3 & th <= th4;
    en(select) = en3 + (en4 - en3) .* rsw(th(select), th3, th4, a34);

    select = th >= th4 & th <= th5;
    en(select) = en5 + (en4 - en5) .* rsw(th(select), th5, th4, a45);

    select = th >= th5 & th <= th6;
    en(select) = en5 + (en6 - en5) .* rsw(th(select), th5, th6, a56);

    select = th >= th6 & th <= th7;
    en(select) = en7 + (en6 - en7) .* rsw(th(select), th7, th6, a67);

end

function en = twists110(th, pars)
    % en = twists110(th,pars)
    %
    % See comments on set110.

    th1 = pars(22); % 110 twist peak position

    en1 = pars(23); % 110 twist energy peak value
    en2 = pars(24); % Sigma3 energy (110 twist, so not a coherent twin)
    en3 = pars(25); % energy at the symmetry point

    a01 = 0.5;
    a12 = 0.5;
    a23 = 0.5;
    %
    th2 = acos(1/3); % Sigma3
    th3 = pi / 2; % 110 90-degree boundary is semi-special, although not a CSL

    perio = pi; % the twist period
    %
    th = mod(abs(th), perio); % rotation symmetry
    th(th > perio / 2) = perio - th(th > perio / 2);

    en = zeros(size(th));
    %
    select = th <= th1;
    en(select) = en1 * rsw(th(select), 0, th1, a01);

    select = th > th1 & th <= th2;
    en(select) = en2 + (en1 - en2) * rsw(th(select), th2, th1, a12);

    select = th > th2;
    en(select) = en3 + (en2 - en3) * rsw(th(select), th3, th2, a23);

end

function en = set111(geom111, pars)
    % en = set111(geom111,pars)
    %
    % Dimensionless contribution to energy from <111> rotations
    % Very similar to set100; see comments therein for general information.
    % Comments in this file will be limited to 111-specific information.

    a = pars(35); % linear part of 111 tilt/twist interpolation
    b = a - 1; % Ensures correct value at x = 1.

    ksi = geom111(2, :);
    eta = geom111(3, :);
    phi = geom111(4, :);

    entwist = twists111(ksi, pars);
    entilt = atgbs111(eta, ksi, pars);
    x = phi / (pi / 2);

    % This one fit well enough with a simple one-parameter parabola that the
    % more complicated power laws in the other sets weren't needed.
    en = entwist + (entilt - entwist) .* (a * x - b * x.^2);
end

function en = twists111(theta, pars)
    % en = twists111(theta,pars)
    %
    % See comments on set111

    thd = pars(37); % 111 twist peak position

    enm = pars(38); % 111 twist energy at the peak
    en2 = pars(28); % Coherent sigma3 twin shows up in two distinct places in the code

    a1 = pars(36); % 111 twist rsw shape parameter
    a2 = a1;

    theta(theta > pi / 3) = 2 * pi / 3 - theta(theta > pi / 3);

    select = (theta <= thd);
    en = zeros(size(theta));
    en(select) = enm * rsw(theta(select), 0, thd, a1);
    en(~select) = en2 + (enm - en2) * rsw(theta(~select), pi / 3, thd, a2);
end

function en = atgbs111(eta, ksi, pars)
    % en = atgbs111(eta,ksi,pars)
    %
    % This function is a fit to the energies of all 111-tilt boundaries

    % There's an additional symmetry in 111 atgbs that doesn't exist in 100 or
    % 110 atgbs.  This is because a rotation about [111] equal to half the period
    % (i.e. 60 degrees) is equivalent to a mirror reflection in the (111)
    % plane.  Both are Sigma3 operations.  The same is not true of the
    % 45-degree [100] or the 90-degree [110] rotation.
    % The following two lines account for this extra symmetry.
    ksi(ksi > pi / 3) = 2 * pi / 3 - ksi(ksi > pi / 3);
    eta(eta > pi / 3) = 2 * pi / 3 - eta(eta > pi / 3);

    % Below the following value of ksi, we ignore the eta dependence.  This is
    % because there's little evidence that it actually varies.  Above this
    % value, we interpolate on an rsw function that follows the Sigma3 line,
    % which is also a line of symmetry for the function.
    ksim = pars(39); % 111 atgb ksi break

    enmax = pars(40); % Energy at the peak (ksi == ksim)
    enmin = pars(41); % energy at the minimum (sigma3, eta == 0)
    encnt = pars(42); % energy at the symmetry point (sigma3, eta == pi/3)

    a1 = 0.5;
    a2 = 0.5;

    etascale = pars(43); % eta scaling parameter for 111 atgb rsw function on Sigma3 line
    % This rsw function is unusual in that the change in shape of the
    % function is much better captured by changing the angular scale rather
    % than changing the dimensionless shape factor.

    en = zeros(size(ksi));

    select = (ksi <= ksim);
    en(select) = enmax * rsw(ksi(select), 0, ksim, a1);

    % chi is the shape of the function along the sigma3 line.
    chi = enmin + (encnt - enmin) * rsw(eta(~select), 0, pi / (2 * etascale), 0.5);
    en(~select) = chi + (enmax - chi) .* rsw(ksi(~select), pi / 3, ksim, a2);
end

function en = rsw(theta, theta1, theta2, a)
    % en = rsw(theta,theta1,theta2,a)
    %
    % This function computes the value of Read-Shockley-Wolf function at theta.
    % The RSW function is normalized to be 1.0 at theta2 and 0.0 at theta1.
    %
    % theta             angle at which to compute the function
    % theta1            the starting angle of the interval
    % theta2            the end angle of the interval
    % a                 parameter defining the shape of the RSW function
    %
    dtheta = theta2 - theta1; % Interval of angles where defined
    theta = (theta - theta1) ./ dtheta * pi / 2; % Normalized angle
    % The rest is the RSW function evaluation
    sins = sin(theta);
    xlogx = zeros(size(sins));

    % Cut off at small sins to avoid 0*infinity problem.  The proper limit is 0.
    select = sins >= 0.000001;
    xlogx(select) = sins(select) .* log(sins(select));

    en = sins - a * xlogx;
end

function [dmin, o2minsyms] = GBdist4(o1, o2, pgnum, dtype, wtol, waitbarQ, epsijk, nv)

    arguments
        o1(:, 8) double {mustBeFinite, mustBeReal, mustBeSqrt2Norm}
        o2(:, 8) double {mustBeFinite, mustBeReal, mustBeSqrt2Norm}
        pgnum(1, 1) double {mustBeInteger} = 32 % default == cubic Oh point group
        dtype char {mustBeMember(dtype, {'omega', 'norm'})} = 'norm'
        wtol double = [] %omega tolerance
        waitbarQ logical = false
        epsijk(1, 1) double = 1
        nv.prec(1, 1) double = 12
        nv.tol(1, 1) double = 1e-6
        nv.nNN(1, 1) double = 1 %number of NNs
        nv.IncludeTies(1, 1) {mustBeLogical} = true
        nv.deleteDuplicates(1, 1) {mustBeLogical} = true
    end

    % GBDIST4  modified version of GBdist function by CMU group. Keeps o1 constant.
    %--------------------------------------------------------------------------
    % Inputs:
    %		o1, o2 -	octonions
    %		pgnum - point group number
    %		dtype - distance type ('omega' arc length or euclidean 'norm')
    % Outputs:
    %		dmin -	minimized distance metric
    %		o2minsyms - minimized octonions
    % Usage:
    %		[dmin, o2minsyms] = GBdist4(o1,o2);
    %		[dmin, o2minsyms] = GBdist4(o1,o2,32);
    %		[dmin, o2minsyms] = GBdist4(o1,o2,32,'norm');
    %		[dmin, o2minsyms] = GBdist4(o1,o2,32,'omega');
    %
    % Dependencies:
    %		osymsets.m
    %			--osymset.m
    %				--qmult.m
    %		get_omega.m
    %		zeta_min2.m (naming distinct from 'zeta_min' to prevent conflicts in
    %		GBdist.m)
    % Notes:
    % Author: Sterling Baird
    % Date: 2020-07-27
    %--------------------------------------------------------------------------
    prec = nv.prec; %precision for duplicates
    tol = nv.tol; %tolerance for duplicates
    nNN = nv.nNN; %number of nearest neighbors
    IncludeTies = nv.IncludeTies; %used in knnsearch
    deleteDuplicatesQ = nv.deleteDuplicates;

    if ~isempty(wtol)
        warning('GBdist4.m functionality has been adjusted to not use wtol (2021-04-14). Specify "prec" and "tol" instead.')
    end

    %number of octonion pairs
    npts = size(o2, 1);

    if (size(o1, 1) == 1) && (size(o2, 1) > 1)
        o1 = repmat(o1, npts, 1);
    end

    grainexchangeQ = true;
    doublecoverQ = true;
    uniqueQ = false;
    %get symmetric octonions (SEOs)
    % osets = osymsets(o2,pgnum,struct,grainexchangeQ,doublecoverQ); %out of memory 2020-08-03

    %assign distance fn to handle
    switch dtype
        case 'omega'
            distfn = @(o1, o2) get_omega(o1, o2);
        case 'alen'
            distfn = @(o1, o2) get_alen(o1, o2);
        case 'norm'
            distfn = @(o1, o2) vecnorm(o1 - o2, 2, 2);
    end

    dmin = zeros(npts, 1);
    o2minsyms = cell(1, npts);

    %textwaitbar setup
    lastwarn('')
    [~, warnID] = lastwarn();
    [~] = ver('parallel');

    if ~strcmp(warnID, 'MATLAB:ver:NotFound')

        try
            D = parallel.pool.DataQueue;
            afterEach(D, @nUpdateProgress);
        catch
            warn("Could start parallel pool")
        end

    else
        waitbarQ = false;
    end

    nsets = npts;
    ninterval = 20;
    N = nsets;
    p = 1;
    reverseStr = '';

    if nsets > ninterval
        nreps2 = floor(nsets / ninterval);
        nreps = nreps2;
    else
        nreps2 = 1;
        nreps = nreps2;
    end

    function nUpdateProgress(~)
        percentDone = 100 * p / N;
        msg = sprintf('%3.0f', percentDone); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        p = p + nreps;
    end

    %loop through octonion pairs, could be sped up significantly via batch approach and/or via GPU adaptation (see gpuArray)
    parfor i = 1:npts %parfor compatible
        %text waitbar
        if mod(i, nreps2) == 0

            if waitbarQ
                send(D, i);
            end

        end

        %% setup
        %unpack symmetrically equivalent octonions (SEOs) of single octonion
        oset = osymsets(o2(i, :), pgnum, struct, grainexchangeQ, doublecoverQ, uniqueQ, epsijk);
        o2tmp = oset{1};

        %number of SEOs
        nsets = size(o2tmp, 1);

        %unpack first octonion (held constant)
        o1tmp = o1(i, :);

        %copy octonion
        o1rep = repmat(o1tmp, nsets, 1);

        %unpack quaternions
        qSC = o2tmp(:, 1:4);
        qSD = o2tmp(:, 5:8);

        %% apply U(1) symmetry
        % get minimum zeta & sigma values (zm)
        zm = zeta_min2(o1rep, o2tmp, -epsijk);
        mA = [0 0 1]; %octonion convention that BP normal is [0 0 1] in lab frame
        mArep = repmat(mA, nsets, 1);
        qzm = ax2qu([mArep zm], -epsijk);
        % 	qzm = [cos(zm/2) zeros(nsets,2) sin(zm/2)];
        %     qzm = [cos(zm/2) zeros(nsets,2) -epsijk*sin(zm/2)];

        % get minimized quaternions
        % 	qCz = qmult(qSC,qzm,epsijk);
        % 	qDz = qmult(qSD,qzm,epsijk);
        qCz = qmult(qzm, qSC, epsijk);
        qDz = qmult(qzm, qSD, epsijk);

        %package quaternions
        o2syms = [qCz qDz];

        %% compute distances
        %give the octonions a norm of sqrt(2)
        o1rep = sqrt2norm(o1rep, 'oct');

        % 	o1rep = repelem(o1rep,8,1);

        %compute all distances
        dlist = distfn(o1rep, o2syms); %#ok<PFBNS> %either omega or euclidean norm (see disttype arg)

        %% find minimum distances & octonions
        %get first instance of minimum omega
        dmin(i) = min(dlist);

        idx = knnsearch(round(dlist, prec), dmin(i), 'IncludeTies', IncludeTies, 'K', nNN);

        if IncludeTies
            minIDs = horzcat(idx{:});
        else
            minIDs = idx;
        end

        %find logical indices of all minimum omegas
        % 	minIDs = ismembertol(dlist,dmin(i),wtol,'DataScale',1); %loosened tol for min omegas, 2020-07-28

        %find corresponding symmetrized octonions (with duplicates)
        o2minsymsTmp = o2syms(minIDs, :);

        %delete duplicate rows (low tol OK b.c. matching 8 numbers)
        if deleteDuplicatesQ
            [~, minIDs] = uniquetol(round(o2minsymsTmp, prec), tol, 'ByRows', true, 'DataScale', 1); %generally repeats will be the same within ~12 sig figs
            o2minsyms{i} = o2minsymsTmp(minIDs, :);
        else
            o2minsyms{i} = o2minsymsTmp; %2021-04-14
        end

    end

end %GBdist4.m

function alen = get_alen(pts, pts2)

    arguments
        pts double
        pts2 double
    end

    % GET_ALEN  get arclength of points on a sphere (acos of dot product, note that this is half the octonion distance)
    pts = normr(pts);
    pts2 = normr(pts2);
    npts = size(pts2, 1);

    if size(pts, 1) == 1
        pts = repmat(pts, npts, 1);
    end

    alen = acos(dot(pts, pts2, 2));
    thr = 1e-6;
    maxim = max(abs(imag(alen)));
    assert(maxim <= thr, ['max imaginary component alen > ' num2str(thr) ', == ' num2str(maxim)])
    alen = real(alen);
end

function q = get_cubo(n, method, sidelength, printQ)

    arguments
        n {mustBeNonNegIntegerOrEmpty} = 1
        method char {mustBeMember(method, {'random', 'uniform'})} = 'random'
        sidelength {mustBeNonNegIntegerOrEmpty} = double.empty
        printQ(1, 1) logical = false
    end

    % GET_CUBO  get n quaternions from randomly or uniformly sampled cubochoric points
    %--------------------------------------------------------------------------
    % Inputs:
    %		n - # of quaternions to output (re-calculated if using 'uniform' method
    %		and sidelength is specified
    %
    %		method - sampling method, 'random', 'uniform'
    %
    %		sidelength - # of points along edge of cube used in cubochoric
    %		sampling (automatically calculated if n is given and sidelength is
    %		not specified)
    %
    % Outputs:
    %		q - rows of quaternions
    %
    % Dependencies: cu2qu.m (De Graef CMU group, see GB Octonion code)
    %
    % References:
    % [1] Singh, S., & De Graef, M. (2016). Orientation sampling for
    % dictionary-based diffraction pattern indexing methods. Modelling and
    % Simulation in Materials Science and Engineering, 24(8).
    % https://doi.org/10.1088/0965-0393/24/8/085013
    %
    % Author: Sterling Baird
    %
    % Date: 2020-07-25
    %--------------------------------------------------------------------------
    if strcmp(method, 'uniform') && isempty(sidelength)
        sidelength = ceil(n^(1/3)); % auto-calculate sidelength
        n = sidelength^3;
    elseif isempty(n)
        n = sidelength^3;
    end

    acube = 0.5 * pi^(2/3);

    %define grid
    switch method
        case 'random'
            G = rand(n, 3); %grid points

        case 'uniform'
            x1 = linspace(0, 1, sidelength);
            [X, Y, Z] = ndgrid(x1, x1, x1);
            G = [X(:) Y(:) Z(:)]; %grid points
    end

    aa = acube * (2 * G - 1); %center grid about [0,0,0] and scale grid

    %convert to quaternion
    q = zeros(n, 4);

    for i = 1:n
        q(i, :) = cu2qu(aa(i, :), printQ); %vectorization possible
    end

end

%fortran code from https://github.com/EMsoft-org/EMsoft/GBs/EMBGO.f90
% acube = 0.5D0 * cPi**0.666666666
%aa = acube * (/ 2*rng_uniform(seed)-1.D0, 2.D0*rng_uniform(seed)-1.D0, 2.D0*rng_uniform(seed)-1.D0 /)
% qa = cu2qu(aa)

%----------------------CUSTOM VALIDATION FUNCTIONS-------------------------
function mustBeNonNegIntegerOrEmpty(arg)
    errmsg = 'must be non-neg integer or empty []';

    if ~isempty(arg)

        if floor(arg) ~= arg
            error(errmsg)
        elseif arg < 0
            error(errmsg)
        end

    end

end

function errout = get_errmetrics(ypred, ytrue, type, NV)

    arguments
        ypred double {mustBeReal}
        ytrue double {mustBeReal} %or measured
        type char {mustBeMember(type, {'e', 'ae', 'mae', 'se', 'rmse', 'all'})} = 'all'
        NV.dispQ(1, 1) logical = false
    end

    % GET_ERRMETRICS  get various error metrics for measured data relative to
    % true data, and specify types based on the lowercase symbol. e.g. 'rmse',
    % or 'all' to output a struct of error metrics. e = predicted - measured
    %--------------------------------------------------------------------------
    % Inputs:
    %  predicted - scalar, vector, or matrix of predicted values
    %
    %  measured - scalar, vector, or matrix of measured (i.e. true or lookup)
    %  values, matching size of predicted
    %
    %  type - type of error metric to output
    %
    % Outputs:
    %  errout - the error metric (or a structure of error metrics) specified by
    %  type, where e, ae, and se, are column vectors,
    %
    % Usage:
    %  errmetrics = get_errmetrics(datatrue,dataout); %outputs struct
    %
    %  rmse = get_errmetrics(datatrue,dataout,'rmse'); %outputs rmse value
    %
    % Dependencies:
    %  var_names.m
    %
    % Notes:
    %  See https://en.wikipedia.org/wiki/Root-mean-square_deviation
    %
    % Author(s): Sterling Baird
    %
    % Date: 2020-09-17
    %--------------------------------------------------------------------------

    %additional argument validation
    szpred = size(ypred);
    szmeas = size(ytrue);
    assert(all(szpred == szmeas), ['predicted size: ' num2str(szpred), ', measured size: ' num2str(szmeas)])

    %error metrics
    e = ypred - ytrue; %error
    ae = abs(e); %absolute error
    mae = mean(ae, 'all'); %mean absolute error
    se = e.^2; %square error
    mse = mean(se, 'all'); %mean square error
    rmse = sqrt(mse); %root mean square error

    if NV.dispQ
        disp(['rmse = ' num2str(rmse) ', mae = ' num2str(mae)])
    end

    %compile into struct
    errmetrics = var_names(ypred, ytrue, e, ae, mae, se, mse, rmse);

    %assign output error value(s)
    switch type
        case 'all'
            errout = errmetrics;
        otherwise
            errout = errmetrics.(type); %dynamic field reference
    end

end

function [datainterp, databary, facetprops, facetIDs, barypars] = ...
        get_interp(mesh, data, intfacetIDs, barytype, barytol, NV)

    arguments
        mesh struct {mustContainFields(mesh, {'pts', 'ppts', 'props', 'sphK'})}
        data struct {mustContainFields(data, {'pts', 'ppts', 'props'})}
        intfacetIDs = []
        barytype char {mustBeMember(barytype, {'planar', 'spherical'})} = 'planar'
        barytol double = double.empty
        NV.saveQ logical = false
        NV.savename char = 'temp.mat'
    end

    % GET_INTERP  Interpolate query point values based on spherical or planar barycentric coords in a mesh
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-27
    %
    % Inputs:
    %	mesh - struct containing octonions (pts), down-projected octonions
    %	(ppts), property values (props), a triangulation (sphK), and a filename
    %	(fname) for the "predictors" and "predictor responses"
    %
    %	data - struct containing octonions (pts), down-projected octonions
    %	(ppts), property values which can be all NaNs (props), and a (fname)
    %	for the query points and interpolated values
    %
    %   intfacetIDs - intersecting facet IDs of data.ppts relative to mesh.ppts
    %   and mesh.sphK
    %
    %   barytype - type of barycentric interpolation, 'spherical' or 'planar'
    %
    %   barytol - tolerance for barycentric coordinates, defaults to 0.2 for
    %   spherical and 1e-6 for planar. *
    %
    % Outputs:
    %	datainterp - interpolated property values at the query points
    %
    %   databary - barycentric coordinates of query points relative to mesh
    %
    %   savename - catenatation of mesh.fname and data.fname where the
    %   workspace is saved
    %
    %   varargout{1} - nndistList, NN euclidean distances
    %
    %   varargout{2} - nonintDists, NN euclidean distances for non-intersecting
    %   points
    %
    % Usage:
    %   [datainterp,databary,savename,varargout] =
    %   get_interp(mesh,data,intfacetIDs,barytype,barytol);
    %
    % Dependencies:
    %   get_omega.m
    %
    %   projray2hypersphere.m
    %
    %   sphbary.m
    %    -projfacet2hyperplane.m
    %     --projray2hyperplane.m
    %
    %   sqrt2norm.m
    %
    % Notes:
    %	consider adding option for 'bary' vs. 'knn' interpolation
    %
    %   *0.2 was chosen as the spherical barytol default because it decreased
    %   the interpolation error relative to planar, whereas lower tolerances
    %   produced more non-intersections.
    %--------------------------------------------------------------------------
    if isempty(barytol)
        %assign default barytol
        switch barytype
            case 'spherical'
                barytol = 0.2;
            case 'planar'
                barytol = 1e-6;
        end

    end

    if isempty(intfacetIDs)
        intfacetIDs = {[]};
    end

    %unpack
    meshpts = mesh.ppts;
    datapts = data.ppts;

    %check dimensions
    if size(mesh.ppts, 2) ~= size(data.ppts, 2)
        errmsg = ['mesh.ppts and data.ppts dims should be equal, but mesh dim == ' ...
                int2str(size(mesh.ppts, 2)) ' and data dim == ' int2str(size(data.ppts, 2))];
        error(errmsg)
    end

    %save name
    % savename = ['mesh_' mesh.fname(1:end-4) '_data_' data.fname];

    %nearest neighbor list
    nnList = dsearchn(meshpts, datapts);
    nndistList = get_omega(sqrt2norm(mesh.pts(nnList, :)), sqrt2norm(data.pts));

    %initialize
    [databary, facetIDs, facetprops] = deal(nan(size(datapts)));
    [datainterp, nonintDists] = deal(nan(size(datapts, 1), 1));
    [nnID, ilist] = deal([]);

    meannorm = mean(vecnorm(mesh.ppts, 2, 2));

    ndatapts = size(datapts, 1);
    disp('loop through datapoints')

    for i = 1:ndatapts
        datapt = datapts(i, :); %use down-projected data (and mesh)
        baryOK = false; %initialize

        if ~isempty(intfacetIDs{i})
            %% setup
            intfacetID = intfacetIDs{i}(1); %take only the first intersecting facet? Average values? Use oSLERP instead?
            vtxIDs = mesh.sphK(intfacetID, :);
            facet = meshpts(vtxIDs, :); %vertices of facet
            facetIDs(i, :) = vtxIDs;
            facetprops(i, :) = mesh.props(vtxIDs).'; %properties of vertices of facet

            %% barycentric coordinates
            switch barytype
                case 'spherical'
                    databary(i, :) = sphbary(datapt, facet); %need to save for inference input
                    nonNegQ = all(databary(i, :) >= -barytol);
                    greaterThanOneQ = sum(databary(i, :)) >= meannorm - barytol;
                    numcheck = all(~isnan(databary(i, :)) & ~isinf(databary(i, :)));
                    baryOK = nonNegQ && greaterThanOneQ && numcheck;

                case 'planar'
                    [~, ~, databaryTemp, ~, ~] = projray2hypersphere(facet, 1:7, datapt, barytol, true);

                    if ~isempty(databaryTemp)
                        databary(i, :) = databaryTemp;
                        nonNegQ = all(databary(i, :) >= -barytol);
                        equalToOneQ = abs(sum(databary(i, :)) - 1) < barytol;
                        numcheck = all(~isnan(databary(i, :)) & ~isinf(databary(i, :)));
                        baryOK = nonNegQ && equalToOneQ && numcheck;
                    end

            end

            if baryOK
                % interpolate using bary coords
                datainterp(i) = dot(databary(i, :), facetprops(i, :));
            else
                disp([num2str(databary(i, :), 2) ' ... sum == ' num2str(sum(databary(i, :)), 10)]);
            end

        end

        if ~baryOK
            disp(['i == ' int2str(i) ...
                    '; no valid intersection, taking NN with dist = ' num2str(nndistList(i))])
            nonintDists(i) = nndistList(i);

            nndistList(i) = NaN; %to distinguish interp vs. NN distances in plotting
            nnID = [nnID nnList(i)]; %#ok<AGROW> %nearest neighbor indices
            ilist = [ilist i]; %#ok<AGROW> % possible to separate out making baryOK a logical array & using 2 for loops
            datainterp(i) = mesh.props(nnList(i)); %assign NN value
        end

    end

    %% Error Metrics
    ids = setdiff(1:ndatapts, ilist); %ids = ~isnan(datainterp);

    interp_errmetrics = get_errmetrics(datainterp(ids), data.props(ids));
    nn_errmetrics = get_errmetrics(mesh.props(nnID), data.props(ilist));

    errmetrics = get_errmetrics( ...
        [datainterp(ids); mesh.props(nnID)], ...
        [data.props(ids); data.props(ilist)]);

    allnn_errmetrics = get_errmetrics(mesh.props(nnList), data.props);

    interpRMSE = interp_errmetrics.rmse;
    nnRMSE = nn_errmetrics.rmse;
    totRMSE = errmetrics.rmse;
    allnnRMSE = allnn_errmetrics.rmse;

    nints = length(ids);
    numnonints = length(ilist);
    int_fraction = nints / (nints + numnonints);

    barypars = var_names(errmetrics, interp_errmetrics, nn_errmetrics, allnn_errmetrics, ...
        ilist, ids, nnList, nndistList, nonintDists, nints, numnonints, int_fraction);

    % fpath = fullfile('data',savename);
    % save(fpath,'-v7.3')

    disp(' ')
    disp(['# non-intersections: ' int2str(sum(~isnan((nnID)))) '/' int2str(ndatapts)])
    disp(' ')
    disp(['RMSE (J/m^2): interp == ' num2str(interpRMSE, '%3.4f') ', NN == ' num2str(nnRMSE, '%3.4f')])
    disp(' ')
    disp(['total RMSE: ' num2str(totRMSE, '%3.4f') ', all NN RMSE comparison: ' num2str(allnnRMSE, '%3.4f')])

    %% Saving
    if NV.saveQ
        save(NV.savename)
    end

end

function [nnX, D, mu, sigma, idxtmp] = get_knn(X, dtype, K, NV)

    arguments
        X double
        dtype char {mustBeMember(dtype, {'omega', 'norm', 'alen'})} = 'omega'
        K(1, 1) double {mustBePositive, mustBeInteger} = 1
        NV.dispQ(1, 1) logical = false
        NV.ID = []
        NV.Y = []
        NV.skipself(1, 1) logical = true
    end

    % GET_KNN  k-nearest neighbor points, distances, mean, and std
    %--------------------------------------------------------------------------
    % Author(s): Sterling Baird
    %
    % Date: 2020-07-27
    %
    % Inputs:
    %  pts - rows of points (euclidean)
    %
    %  dtype - distance type
    %
    % Outputs:
    %  creates a histogram figure
    %
    % Usage:
    %	pts = sqrt2norm(normr(rand(1000,8))); %octonion case (omega), input must have norm == sqrt(2)
    %  nnhist(pts)
    %
    %  pts = nnhist(rand(1000,5)); % euclidean distance
    %  nnhist(pts)
    %
    % Dependencies:
    %  get_omega.m
    %
    % Notes:
    %  *
    %--------------------------------------------------------------------------
    %get nearest neighbor IDs and euclidean distances
    ID = NV.ID;
    Y = NV.Y;
    skipselfQ = NV.skipself;
    assert(isempty(ID) || isempty(Y), 'ID and Y should not be supplied simultaneously')

    if ~isempty(ID)
        Y = X(ID, :);
    elseif isempty(Y)
        Y = X; %check NNs within set of points (otherwise check NN against specific pts in NV.Y)
    end

    if skipselfQ
        [idxtmp, Dtmp] = knnsearch(X, Y, 'K', K + 1);
    else
        [idxtmp, Dtmp] = knnsearch(X, Y, 'K', K);
    end

    %remove "self" column
    if skipselfQ
        idxtmp = idxtmp(:, 2:end);
        Dtmp = Dtmp(:, 2:end);
    end

    %initialize
    [mu, sigma] = deal(zeros(K, 1));
    D = cell(K, 1);

    for k = 1:K
        idx = idxtmp(:, k);

        %nearest neighbor pts
        nnX = X(idx, :);

        %get distances for plotting
        switch dtype
            case 'omega'
                Drad = get_omega(nnX, Y);
                D{k} = rad2deg(Drad);
            case 'norm'
                D{k} = Dtmp(:, k);
            case 'alen' %general arc length formula
                assert(norm(X(1, :)) - 1 < 1e-6, 'norm(pts,2) must be 1 (within tol)')
                Drad = real(acos(dot(nnX, Y, 2)));
                D{k} = rad2deg(Drad);
        end

        mu(k) = mean(D{k});
        sigma(k) = std(D{k});

        if NV.dispQ
            disp(['nn: ' int2str(k) ', mu = ' num2str(mu(k)) ', sigma = ' num2str(sigma(k))])
        end

    end

end

function [octvtx, oref, fiveref, ids] = get_octpairs(pts, epsijk, nv)

    arguments
        pts(:, 8) double {mustBeSqrt2Norm}
        epsijk(1, 1) double = 1
        nv.o2addQ(1, 1) logical = false
        nv.pgnum(1, 1) double = 32
        nv.wtol double = []
        nv.fiveref = []
        nv.oref(1, 8) double = get_ocubo(1, 'random', [], 10)
        nv.dispQ = []
        nv.nNN(1, 1) double = 1 %number of NNs
        nv.IncludeTies(1, 1) {mustBeLogical} = true
    end

    % GET_OCTPAIRS  Get a set of octonions that are symmetrized with respect to a fixed reference GB (default rng seed == 10)
    % Author: Sterling Baird
    %
    % Date: 2020-07-27
    %
    % Inputs:
    %  pts - rows of octonions
    %
    %  o2addQ -	logical, whether to add o2 (in this case, point 'O' with
    %		nA = [0 0 1]) to "five", and if not, then get rid
    %       of pts(1,:) which corresponds to the same.
    %
    % Outputs:
    %  octvtx - rows of octonions that form a mesh
    %
    %  usv - struct to use with proj_down.m and proj_up.m
    %
    %  five - struct containing misorientation quaternions (q), BP normals (nA,
    %  grain A reference frame), rodrigues vectors (d), and misorientation
    %  fundamental zone feature type (geometry)

    % Dependencies:
    %  misFZfeatures.mat
    %  GBdist4.m
    %  mustBeSqrt2Norm.m (argument validation function)
    %--------------------------------------------------------------------------
    dispQ = nv.dispQ;
    nNN = nv.nNN;
    IncludeTies = nv.IncludeTies;

    if isempty(dispQ)

        if size(pts, 1) <= 1000
            dispQ = false;
        else
            dispQ = true;
        end

    end

    %% Unpack 5DOF reference (empty is OK)
    fiveref = nv.fiveref;

    %% get reference octonion
    if isempty(fiveref)
        oref = nv.oref;
    else
        oref = five2oct(fiveref, epsijk);
    end

    %% get minimized distance octonions relative to oct pairs
    if dispQ
        disp('get_octpairs ')
    end

    npts = size(pts, 1);
    % orefrep = repmat(oref,npts,1);
    [dmin, octvtx] = GBdist4(oref, pts, nv.pgnum, 'norm', nv.wtol, dispQ, epsijk, 'IncludeTies', IncludeTies, 'nNN', nNN);

    len = cellfun(@(x) size(x, 1), octvtx);
    ids = arrayfun(@(x) repelem(x, len(x)), 1:npts, 'UniformOutput', false);
    ids = [ids{:}];

    %catenate
    octvtx = vertcat(octvtx{:});

    if nv.o2addQ
        %add reference octonion
        octvtx = [oref; octvtx];
    end

end

function o = get_ocubo(n, method, sidelength, seed)

    arguments
        n {mustBeNonNegIntegerOrEmpty} = 1
        method char {mustBeMember(method, {'random', 'uniform'})} = 'random'
        sidelength {mustBeNonNegIntegerOrEmpty} = []
        seed = []
    end

    % GET_OCUBO  get octonions formed by pairs of quaternions from randomly or uniformly sampled cubochoric points.
    %  In general, for random, no two quaternions will be the same.
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-25
    %
    % Inputs:
    %		n - # of octonions to output (re-calculated if using 'uniform' method
    %		and sidelength is specified
    %
    %		method - sampling method, 'random', 'uniform'
    %
    %		sidelength - # of points along edge of cube used in cubochoric
    %		sampling (automatically calculated if n is given and sidelength is
    %		not specified)
    %
    % Outputs:
    %
    %		o - list of octonions
    %
    % Usage:
    %		o = get_ocubo(); %generate a single octonion formed by two
    %		quaternions sampled randomly from cubochoric space
    %
    %		o = get_ocubo(5) %generate 5 octonions from pairs of quaternions
    %		randomly sampled from cubochoric space
    %
    %		o = get_ocubo(5,'random') %generate 5 octonions from pairs of
    %		quaternions randomly sampled from cubochoric space
    %
    %		o = get_ocubo(100,'uniform') %generate 100 octonions randomly sampled
    %		from list of pairs of quaternions generated via uniform cubochoric
    %		sampling with automatically calculated sidelength (ceil(100^(1/3))
    %
    %		o = get_ocubo([],'uniform',5) %generate all combinations of
    %		quaternion pairs (i.e. octonions) using 5^3 == 125 uniformly sampled
    %		quaternions (15625 octonions)
    %
    %       o = get_ocubo(100,'random',[],10) %generate 100 random octonions
    %       using a random number generator seed of 10
    %
    % Dependencies:
    %		allcomb.m (optional if nboQ == false)
    %
    %		ocubo.m
    %			--cu2qu.m (De Graef group)
    %
    % Note: specifying 'random' and a non-empty sidelength will error, as these
    % are two contradictory options.
    %
    %--------------------------------------------------------------------------
    %set random number generator (only if custom seed is specified)
    if ~isempty(seed)
        startseed = rng; %to be used at end to set rng to initial state (i.e. before calling this function)
        rng(seed)
    end

    % argument validation (cont.)
    if strcmp(method, 'random') && ~isempty(sidelength)
        error('sidelength should not be specified for random sampling')
    end

    %--------------------------------------------------------------------------
    %setup

    if strcmp(method, 'uniform') && isempty(sidelength)
        sidelength = ceil(n^(1/3)); % auto-calculate sidelength
        nq = [];
    else

        if isempty(sidelength)
            nq = n;
        else
            nq = sidelength^3; %auto-calculate # of quaternions
        end

    end

    %--------------------------------------------------------------------------

    switch method
        case 'random'
            % get 2*n cubochoric points
            q = get_cubo(2 * nq, method, sidelength);

            %unpack
            qA = q(1:n, :);
            qB = q(n + 1:2 * n, :);

        case 'uniform'
            %get n cubochoric points
            q = get_cubo(nq, method, sidelength);

            %convert to cell array of quaternions
            q = num2cell(q, 2);

            %form pairs
            nboQ = true; %whether to include no-boundary octonions

            if nboQ
                qpairs = allcomb(q, q);
            else
                qpairs = nchoosek(q, 2);
            end

            % unpack
            qA = vertcat(qpairs{:, 1});
            qB = vertcat(qpairs{:, 2});

            if ~isempty(n)

                if (n < length(q)) && ~isempty(sidelength)
                    % get a random list of the quaternions
                    randlist = randi(size(qA, 1), n, 1);
                    qA = qA(randlist, :);
                    qB = qB(randlist, :);
                end

            end

    end

    %catenate
    o = [qA qB];

    %get unique list
    [~, ia] = uniquetol(round(o, 12), 'ByRows', true);
    o = o(ia, :);

    if ~isempty(seed)
        %reset rng back to what it was before calling the function
        rng(startseed);
    end

end %get_ocubo.m

function omega = get_omega(o1, o2)

    arguments
        o1 (:, 8) double {mustBeReal, mustBeFinite}
        o2 (:, 8) double {mustBeReal, mustBeFinite}
    end

    % GET_OMEGA  calculate octonion distance
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date:
    %
    % Description:
    %
    % Inputs: o1, o2 form an octonion pair, each can be rows of octonions
    %
    % Outputs:
    %
    % Dependencies:
    %
    %--------------------------------------------------------------------------

    % enables use as custom distance function to e.g. pdist2
    if size(o1, 1) == 1
        npts = size(o2, 1);
        o1 = repmat(o1, npts, 1);
    end

    if size(o2, 1) == 1
        npts = size(o1, 1);
        o2 = repmat(o2, npts, 1);
    end

    qA = normr(o1(:, 1:4));
    qB = normr(o1(:, 5:8));
    qC = normr(o2(:, 1:4));
    qD = normr(o2(:, 5:8));

    dot1 = dot(qA, qC, 2);
    dot2 = dot(qB, qD, 2);

    omega = real(2 * acos(abs(dot1 + dot2) / 2)); %added real() 2020-08-03

end

function [Spairs, nsympairs] = get_sympairs(pgnum, dispnameQ)

    arguments
        pgnum(1, 1) int8 {mustBeInteger} = 32 %default == cubic
        dispnameQ(1, 1) logical = false
    end

    % GET_SYMPAIRS  get all combinations (pairs) of operators for a point group
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-27
    %
    % Inputs:
    %		pgnum - point group number, 32 == Oh (cubic)
    %
    %		dispnameQ - whether or not to display the name of the point group
    %
    % Outputs:
    %		Spairs - combinations of symmetry operators, including repeats
    %
    % Usage:
    %		Spairs = get_sympairs(32); % combinations of cubic symmetry operators
    %
    % Dependencies:
    %		Pgsymops.mat, PGnames.mat (optional, default is not required)
    %
    %--------------------------------------------------------------------------

    %load operators
    % symops = load('PGsymops.mat');
    Q = {[
        1 0 0 0

        ], [

        1 0 0 0

        ], [

        1 0 0 0
        0 0 1 0

        ], [

        1 0 0 0

        ], [

        1 0 0 0
        0 0 1 0

        ], [

        1 0 0 0
        0 0 1 0

        ], [

        1 0 0 0
        0 0 1 0

        ], [

        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1

        ], [

        1.0000 0 0 0
        0 0 0 1.0000
        0.7071 0 0 0.7071
        0.7071 0 0 -0.7071

        ], [

        1.0000 0 0 0
        0 0 0 1.0000
        0.7071 0 0 0.7071
        0.7071 0 0 -0.7071

        ], [

        1.0000 0 0 0
        0 0 0 1.0000
        0.7071 0 0 0.7071
        0.7071 0 0 -0.7071

        ], [

        1.0000 0 0 0
        0 0 0 1.0000
        0.7071 0 0 0.7071
        0.7071 0 0 -0.7071
        0 1.0000 0 0
        0 0 1.0000 0
        0 0.7071 0.7071 0
        0 -0.7071 0.7071 0

        ], [

        1.0000 0 0 0
        0 0 0 1.0000
        0.7071 0 0 0.7071
        0.7071 0 0 -0.7071

        ], [

        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1

        ], [

        1.0000 0 0 0
        0 0 0 1.0000
        0.7071 0 0 0.7071
        0.7071 0 0 -0.7071
        0 1.0000 0 0
        0 0 1.0000 0
        0 0.7071 0.7071 0
        0 -0.7071 0.7071 0

        ], [

        1.0000 0 0 0
        0.5000 0 0 0.8660
        -0.5000 0 0 0.8660

        ], [

        1.0000 0 0 0
        0.5000 0 0 0.8660
        -0.5000 0 0 0.8660

        ], [

        1.0000 0 0 0
        0.5000 0 0 0.8660
        -0.5000 0 0 0.8660
        0 1.0000 0 0
        0 0.5000 0.8660 0
        0 -0.5000 0.8660 0

        ], [

        1.0000 0 0 0
        0.5000 0 0 0.8660
        -0.5000 0 0 0.8660

        ], [

        1.0000 0 0 0
        0.5000 0 0 0.8660
        -0.5000 0 0 0.8660
        0 1.0000 0 0
        0 0.5000 0.8660 0
        0 -0.5000 0.8660 0

        ], [

        1.0000 0 0 0
        0.8660 0 0 0.5000
        0.5000 0 0 0.8660
        0 0 0 1.0000
        -0.5000 0 0 0.8660
        -0.8660 0 0 0.5000

        ], [

        1.0000 0 0 0
        0.8660 0 0 0.5000
        0.5000 0 0 0.8660
        0 0 0 1.0000
        -0.5000 0 0 0.8660
        -0.8660 0 0 0.5000

        ], [

        1.0000 0 0 0
        0.8660 0 0 0.5000
        0.5000 0 0 0.8660
        0 0 0 1.0000
        -0.5000 0 0 0.8660
        -0.8660 0 0 0.5000

        ], [

        1.0000 0 0 0
        0.8660 0 0 0.5000
        0.5000 0 0 0.8660
        0 0 0 1.0000
        -0.5000 0 0 0.8660
        -0.8660 0 0 0.5000
        0 1.0000 0 0
        0 0.8660 0.5000 0
        0 0.5000 0.8660 0
        0 0 1.0000 0
        0 -0.5000 0.8660 0
        0 -0.8660 0.5000 0

        ], [

        1.0000 0 0 0
        0.8660 0 0 0.5000
        0.5000 0 0 0.8660
        0 0 0 1.0000
        -0.5000 0 0 0.8660
        -0.8660 0 0 0.5000

        ], [

        1.0000 0 0 0
        0.8660 0 0 0.5000
        0.5000 0 0 0.8660
        0 0 0 1.0000
        -0.5000 0 0 0.8660
        -0.8660 0 0 0.5000

        ], [

        1.0000 0 0 0
        0.8660 0 0 0.5000
        0.5000 0 0 0.8660
        0 0 0 1.0000
        -0.5000 0 0 0.8660
        -0.8660 0 0 0.5000
        0 1.0000 0 0
        0 0.8660 0.5000 0
        0 0.5000 0.8660 0
        0 0 1.0000 0
        0 -0.5000 0.8660 0
        0 -0.8660 0.5000 0

        ], [

        1.0000 0 0 0
        0 1.0000 0 0
        0 0 1.0000 0
        0 0 0 1.0000
        0.5000 0.5000 0.5000 0.5000
        0.5000 -0.5000 -0.5000 -0.5000
        0.5000 0.5000 -0.5000 0.5000
        0.5000 -0.5000 0.5000 -0.5000
        0.5000 -0.5000 0.5000 0.5000
        0.5000 0.5000 -0.5000 -0.5000
        0.5000 -0.5000 -0.5000 0.5000
        0.5000 0.5000 0.5000 -0.5000

        ], [

        1.0000 0 0 0
        0 1.0000 0 0
        0 0 1.0000 0
        0 0 0 1.0000
        0.5000 0.5000 0.5000 0.5000
        0.5000 -0.5000 -0.5000 -0.5000
        0.5000 0.5000 -0.5000 0.5000
        0.5000 -0.5000 0.5000 -0.5000
        0.5000 -0.5000 0.5000 0.5000
        0.5000 0.5000 -0.5000 -0.5000
        0.5000 -0.5000 -0.5000 0.5000
        0.5000 0.5000 0.5000 -0.5000

        ], [

        1.0000 0 0 0
        0 1.0000 0 0
        0 0 1.0000 0
        0 0 0 1.0000
        0.7071 0.7071 0 0
        0.7071 0 0.7071 0
        0.7071 0 0 0.7071
        0.7071 -0.7071 0 0
        0.7071 0 -0.7071 0
        0.7071 0 0 -0.7071
        0 0.7071 0.7071 0
        0 -0.7071 0.7071 0
        0 0 0.7071 0.7071
        0 0 -0.7071 0.7071
        0 0.7071 0 0.7071
        0 -0.7071 0 0.7071
        0.5000 0.5000 0.5000 0.5000
        0.5000 -0.5000 -0.5000 -0.5000
        0.5000 0.5000 -0.5000 0.5000
        0.5000 -0.5000 0.5000 -0.5000
        0.5000 -0.5000 0.5000 0.5000
        0.5000 0.5000 -0.5000 -0.5000
        0.5000 -0.5000 -0.5000 0.5000
        0.5000 0.5000 0.5000 -0.5000

        ], [

        1.0000 0 0 0
        0 1.0000 0 0
        0 0 1.0000 0
        0 0 0 1.0000
        0.5000 0.5000 0.5000 0.5000
        0.5000 -0.5000 -0.5000 -0.5000
        0.5000 0.5000 -0.5000 0.5000
        0.5000 -0.5000 0.5000 -0.5000
        0.5000 -0.5000 0.5000 0.5000
        0.5000 0.5000 -0.5000 -0.5000
        0.5000 -0.5000 -0.5000 0.5000
        0.5000 0.5000 0.5000 -0.5000

        ], [

        1.0000 0 0 0
        0 1.0000 0 0
        0 0 1.0000 0
        0 0 0 1.0000
        0.7071 0.7071 0 0
        0.7071 0 0.7071 0
        0.7071 0 0 0.7071
        0.7071 -0.7071 0 0
        0.7071 0 -0.7071 0
        0.7071 0 0 -0.7071
        0 0.7071 0.7071 0
        0 -0.7071 0.7071 0
        0 0 0.7071 0.7071
        0 0 -0.7071 0.7071
        0 0.7071 0 0.7071
        0 -0.7071 0 0.7071
        0.5000 0.5000 0.5000 0.5000
        0.5000 -0.5000 -0.5000 -0.5000
        0.5000 0.5000 -0.5000 0.5000
        0.5000 -0.5000 0.5000 -0.5000
        0.5000 -0.5000 0.5000 0.5000
        0.5000 0.5000 -0.5000 -0.5000
        0.5000 -0.5000 -0.5000 0.5000
        0.5000 0.5000 0.5000 -0.5000
        ]};

    if dispnameQ
        %display point group name
        % symnames = load('PGnames.mat');
        PG_names = {'1: 1', '2: -1', '3: 2', '4: m', '5: 2/m', '6: 222', '7: mm2', '8: mmm', '9: 4', '10: -4', '11: 4/m', '12: 422', '13: 4mm', '14: -42m', '15:4/mmm', '16: 3', '17: -3', '18: 32', '19: 3m', '20: -3m', '21: 6', '22: -6', '23:  6/m', '24:  622', '25: 6mm', '26: -6m2', '27: 6/mmm', '28: 23', '29: m3', '30: 432', '31: -43m', '32: m-3m'};
        % pgname = symnames.PG_names{pgnum};
        pgname = PG_names{pgnum};
        disp(pgname)
    end

    %unpack point group
    % qpt = symops.Q{pgnum};
    qpt = Q{pgnum};

    %get all combinations of symmetry operators
    qpttmp = num2cell(qpt, 2);
    qptpairs = allcomb(qpttmp, qpttmp);

    % unpack symmetry operator combinations
    SAlist = vertcat(qptpairs{:, 1});
    SBlist = vertcat(qptpairs{:, 2});

    % catenate combinations
    Spairs = [SAlist SBlist];

end %get_sympairs

function [yq, W, r, nints, numnonints, int_fraction] = idw(X, Xq, y, r, L)

    arguments
        X double %rows of input points
        Xq double %rows of query points
        y(:, 1) double %property values
        r double = [] %radius
        L(1, 1) double = 2 %Euclidean or 2-norm
    end

    % IDW inverse-distance weighting interpolation
    %referenced Andres Tovar's FEX implementation (idw.m)

    %default radius or user-supplied radius
    if isempty(r)
        [~, ~, mu] = get_knn(X, 'norm', 1);
        r = mu * sqrt(2); %note: since alen = ~0.5*omega, mu in octonion is ~2x this (2*mu*sqrt(2))
    end

    pd = pdist2(X, Xq); %pairwise distance matrix

    %parse pd based on radius
    pd(pd == 0) = eps;
    pd(pd > r) = Inf;

    W = 1 ./ (pd.^L); %weight matrix

    %store as sparse matrix if more than half of elements are zero
    numnonzero = nnz(W);

    if numnonzero < 0.5 * numel(pd)
        W = sparse(W);
    end

    %compute interpolated values
    yq = sum(W .* y) ./ sum(W);

    %assign NN value if no GBs fall into radius for a particular property
    nnIDs = find(isnan(yq));
    nnList = dsearchn(X, Xq(nnIDs, :));
    %assign NN property values
    yq(nnIDs) = y(nnList);

    if issparse(yq)
        yq = full(yq);
    end

    %make into column vector
    yq = yq.';

    %idw parameters
    numnonints = length(nnIDs);
    npredpts = size(Xq, 1);
    nints = npredpts - numnonints;
    int_fraction = nints / (nints + numnonints);
    disp(['numnonints = ' int2str(numnonints)])

end

function [intfacetIDs, dataBary, klist] = intersect_facet(pts, K, datalist, tol, NV)

    arguments
        pts double {mustBeFinite, mustBeReal}
        K uint16 {mustBeFinite, mustBeReal} %maximum of 65535 points (not rows), change to uint32 otherwise
        datalist double {mustBeFinite, mustBeReal}
        tol(1, 1) double {mustBeFinite, mustBeReal} = 1e-6
        NV.maxnormQ(1, 1) logical = true
        NV.inttype char {mustBeMember(NV.inttype, {'planar', 'spherical'})} = 'planar'
        NV.invmethod char {mustBeMember(NV.invmethod, {'mldivide', 'pinv', 'extendedCross'})} = 'mldivide'
        NV.nnMax(1, 1) double {mustBeInteger} = 10 %size(pts,1)
    end

    % INTERSECT_FACET  Find intersection of ray with facet using barycentric coordinates.
    % Project a ray (each row of pts) onto each facet connected to
    % the nearest neighbor of that ray, and compute the barycentric coordinates
    % of the projected datapoint. If all coordinates of a facet are positive,
    % mark that facet as an intersecting facet. If no intersecting is found
    % with the first nearest neighbor search, continue looking at next nearest
    % neighbors until an intersecting facet has been found or all facets have
    % been looped through.
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-06
    %
    % Inputs:
    %		pts			===	rows of points that fall on a unit sphere, must
    %								match indices in K
    %
    %		K				===	convex hull triangulation of pts.
    %
    %		datalist		===	rows of datapoints that represent rays from origin
    %								to the datapoint to be projected onto the facet.
    %
    % Outputs:
    %		intfacetIDs ===	intersecting facet IDs (as in row index of K).
    %								Returns NaN if no intersecting facet is found, even
    %								after looping through all facets instead of just
    %								ones connected to first NN. Facet with
    %								largest norm returned if multiple facets found.
    %
    % Dependencies:
    %		projray2hypersphere.m
    %			numStabBary.m (optional)
    %
    % Notes:
    %
    %--------------------------------------------------------------------------
    maxnormQ = NV.maxnormQ;
    inttype = NV.inttype;
    invmethod = NV.invmethod;
    nnMax = NV.nnMax;

    %% find nearest vertex for each datapoint
    nnList = dsearchn(pts, datalist);

    % nmeshpts = size(pts,1);
    ndatapts = size(datalist, 1);

    klist = zeros(ndatapts, 1);

    %% initalize
    dataProj = cell(1, ndatapts); %projected data
    facetPts = dataProj; %facet points
    dataBary = dataProj; %barycentric data
    subfacetIDs = dataProj; %sub IDs (i.e. IDs from sublist of facets)
    intfacetIDs = dataProj; % facet IDs that can be used to index into K
    t = dataProj;

    waitbarQ = true;
    %textwaitbar setup
    lastwarn('')
    [~, warnID] = lastwarn();
    [~] = ver('parallel');

    if ~strcmp(warnID, 'MATLAB:ver:NotFound')
        D = parallel.pool.DataQueue;
        K = parallel.pool.Constant(K);
        afterEach(D, @nUpdateProgress);
    else
        waitbarQ = false;
    end

    N = ndatapts;
    p = 1;
    reverseStr = '';
    nreps2 = ceil(N / 100);
    nreps = nreps2;

    function nUpdateProgress(~)
        percentDone = 100 * p / N;
        msg = sprintf('%3.0f', percentDone); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        p = p + nreps;
    end

    %% loop through datapts
    parfor i = 1:ndatapts % parfor compatible
        %text waitbar
        if mod(i, nreps2) == 0

            if waitbarQ
                send(D, i);
            end

        end

        %% first NN projection
        data = datalist(i, :);
        nn = nnList(i);

        rownext = []; %initialize (used in while loop)

        %find vertices of facets attached to NN vertex (or use all facets)
        [row, ~] = find(K.Value == nn);
        facetPtIDs = K.Value(row, :);

        %compute projections
        switch inttype
            case 'planar'
                [dataProj{i}, facetPts{i}, dataBary{i}, subfacetIDs{i}, t{i}] = ...
                    projray2hypersphere(pts, facetPtIDs, data, tol, maxnormQ, invmethod);
            case 'spherical'
                [dataBary{i}, subfacetIDs{i}] = sphbary_setup(pts, facetPtIDs, data, tol); %I think this is buggy 2020-07-16
        end

        %% keep using next NNs if facet not found
        ptsTemp = pts; %dummy variable to be able to sift through new NN's
        k = 0;
        oldrow = row;

        while isempty(subfacetIDs{i}) && k < nnMax
            k = k + 1;
            %remove previous NN
            ptsTemp(nn, :) = NaN(1, size(pts, 2));

            %find next NN
            nn = dsearchn(ptsTemp, data);

            %find facets attached to next NN
            [row, ~] = find(K.Value == nn);

            rownext = setdiff(row, oldrow); %this seems problematic for indexing later (2020-07-29)
            oldrow = [row; oldrow];

            if ~isempty(rownext)
                facetPtIDsNext = K.Value(rownext, :);

                %compute projections
                switch inttype
                    case 'planar'
                        [dataProj{i}, facetPts{i}, dataBary{i}, subfacetIDs{i}, t{i}] = ...
                            projray2hypersphere(pts, facetPtIDsNext, data, tol, maxnormQ, invmethod);
                    case 'spherical'
                        [dataBary{i}, subfacetIDs{i}] = sphbary_setup(pts, facetPtIDs, data, tol);
                end

            end

        end

        if k > 0
            row = rownext; %correct for indexing if the while loop was entered into
            %NOTE: not having this was a major source of error (2020-07-29),
            %i.e. only datapoints which did not enter the while loop had the
            %correct output for the intersecting facet ID
        end

        if ~isempty(subfacetIDs{i})
            %convert from facetPtIDs or facetPtIDsNext index to K index
            intfacetIDs{i} = row(subfacetIDs{i});
        else
            intfacetIDs{i} = [];
        end

        klist(i) = k;
    end

end

function mustBeLogical(value)
    % MUSTBELOGICAL  must be a logical value, error otherwise
    %NOTE: do not use "arguments" syntax here since this is a validation fn
    %---------------------HELPER VALIDATION FUNCTION---------------------------
    if ~islogical(value)
        error('value must be of type logical')
    end

    %mustBeA(value,'logical')
end

function mustBeSqrt2Norm(o)
    % MUSTBESQRT2NORM  check that first octonion in list has norm == sqrt(2) and each quaternion has norm == 1
    %NOTE: do not use "arguments" syntax here since this is a validation fn
    %---------------------HELPER VALIDATION FUNCTION---------------------------
    onorm = norm(o(1, :));
    errmsg = ['octonion norm == ' num2str(onorm) ' ~= sqrt(2) == 1.4142'];
    assert(abs(onorm - sqrt(2)) <= 1e-3, errmsg)
end %mustBeSqrt2Norm

function mustContainFields(S, checknames)
    % MUSTCONTAINFIELDS  check fieldnames(S) and make sure every checkname exists
    varnames = fieldnames(S);
    errmsg = ['input needs minimum fields: ' strjoin(checknames) ...
            ' but contains fields: ' strjoin(varnames)];
    assert(all(ismember(checknames, varnames)), errmsg)
end

function nrlist = normr(rlist)
    % NORMR  normalizes vectors row-by-row. Outputs zero vector if zero vector input (shadowed by Computer Vision Toolbox)
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-03
    %
    % Inputs:
    %
    %		rlist		===	rows of vectors to be normalized
    %
    % Outputs:
    %
    %		nrlist	===	normalized rows of vectors
    %
    % Dependencies:
    %
    % Note:	normr.m shadows a built-in function that's part of the Computer
    %			Vision Toolbox
    %
    %--------------------------------------------------------------------------

    %determine size
    [n, d] = size(rlist);

    % shortcut if one row
    if n == 1
        nm = norm(rlist);

        if nm ~= 0
            nrlist = rlist ./ norm(rlist);
        else
            nrlist = zeros(n, d);
        end

    else
        %initialize
        nrlist = zeros(size(rlist));

        %compute norms
        nmlist = vecnorm(rlist, 2, 2);

        %get indices of non-zero elements
        ids = find(nmlist);

        %normalize only rows where norm is non-zero
        nrlist(ids, :) = rlist(ids, :) ./ nmlist(ids);

        %note: when nm(id) == 0, nrlist(id) == zeros(1,d)
    end

end %normr

function lambda = numStabBary(X, xA)
    % NUMSTABBARY  a numerically stable barycentric approach in high dimensions
    % ------NOTES------
    %
    % A numerically stable approach to computing barycentric coordinates by
    % converting to projective or homogeneous coordinates as an intermediate to
    % cartesian. I.e. cartesian -> projective -> barycentric. Uses an n-ary cross
    % product [1] computed with MATLAB's symbolic engine to solve a linear system
    % of equations (instead of using `\` (also called mldivide() or backsolve)
    % operator). Taking this approach should theoretically minimize division
    % operations that commonly cause numerical instability [2].
    %
    % Max is 8D cartesian, add more symbolic vectors to U and I for higher
    % dimensions
    %
    % ------INPUT-------
    %
    % X - row-wise vertices of simplex xA - test point that might be within
    % simplex
    %
    % ---INTERMEDIATE---
    %
    % S - rows of vectors defining simplex, with X(1,:) as origin (i.e.
    % parametric representation)
    %
    % ------OUTPUT------
    %
    % lambda - barycentric coordinates
    %
    %
    % ----4D EXAMPLE----
    %
    % X = [ 1 1 0 1; 1 1 1 0; 0 1 1 1; 1 0 1 1 ];
    %
    % xA = [1 1 1 1];
    %
    % lambda = numStabBary(X,xA);
    %
    % ----REFERENCES----
    %
    % [1] https://en.wikipedia.org/wiki/Cross_product
    %
    % [2] V. Skala, Robust Barycentric Coordinates Computation of the Closest
    % Point to a Hyperplane in E^n, Proc. 2013 Int. Conf. Applies Math. Comput.
    % Methods Eng. (2013) 239�244.

    d = length(xA);

    x1 = X(1, :); %first vertex

    S = zeros(size(X) - [1 0]); %one less row than X
    srows = size(S, 1);

    % convert to parametric representation
    for i = 2:size(X, 1)
        S(i - 1, :) = X(i, :) - x1;
    end

    % construct A matrix as in Ax=b
    A = zeros(srows, srows);

    for i = 1:srows

        for j = i:srows
            A(i, j) = dot(S(i, :), S(j, :)); %creates upper triangular matrix
        end

    end

    A = tril(A.') + triu(A, 1); % use upper triangle to make symmetric

    %disp(['det(A) = ',num2str(det(A))]) %gives sense of stability (how thin the simplex is)

    %construct b vector as in Ax=b
    b = zeros(srows, 1);

    for i = 1:srows
        b(i) = dot(S(i, :), x1 - xA);
    end

    method = 'extendedCross'; %'backslash', 'extendedCross'

    switch method
        case 'extendedCross'
            %convert to extended cross product form (n-ary cross product)
            eta = zeros(srows, srows + 1); %one more column to accomodate "b", one more row of 1's

            for i = 1:srows
                eta(i, :) = [A(i, :) b(i)];
            end

            %compute n-ary cross product (sigma)
            syms i1 i2 i3 i4 i5 i6 i7 i8 i9
            I = [i1 i2 i3 i4 i5 i6 i7 i8 i9]; %unit vectors
            eta = [eta; I(1:d)]; % Wikipedia suggests adding to last row so that direction of othorgonal vector is correct. Conflicts with reference in paper above
            etadet = det(eta); %can't just take determinant - top or bottom row needs to be i,j,k,l unit vectors, not 1,1,1,1

            if etadet == 0
                %disp('etadet == 0');
                lambda = repelem(-Inf, d);
                return
            end

            [C, T] = coeffs(etadet, I(1:d)); % k x 1 vector orthogonal to k-1 vectors
            Cvars = ismember(I(1:d), T);
            sigma(Cvars) = C;
            sigma(~Cvars) = 0;
            %compute parametric coordinates of point in Euclidean space

            U = double(sigma(1:end - 1) / sigma(end));
        case 'backslash'
            U = A \ b;
    end

    %compute barycentric coordinates
    lambda = zeros(1, d);
    lambda(1) = 1 - sum(U);

    try
        lambda(2:end) = U;
    catch
        lambda(2:end) = U;
    end

end %numStabBary

function symocts = osymset(qA, qB, Spairs, grainexchangeQ, doublecoverQ, uniqueQ, epsijk)

    arguments
        qA(1, 4) double {mustBeNumeric, mustBeFinite}
        qB(1, 4) double {mustBeNumeric, mustBeFinite}
        Spairs(:, 8) double {mustBeNumeric} = get_sympairs(32, false) %default to cubic Oh symmetry
        grainexchangeQ(1, 1) logical {mustBeLogical} = false
        doublecoverQ(1, 1) logical {mustBeLogical} = false
        uniqueQ(1, 1) logical {mustBeLogical} = false
        epsijk(1, 1) double {mustBeInteger} = 1
    end

    % OSYMSET  get symmetrically equivalent octonions
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-27
    %
    % Inputs:
    %		(qA, qB) - quaternions
    %
    %		Slist - list of pairs of symmetry operators to be applied to qA and qB
    %
    % Outputs: rows of symmetrically equivalent octonions
    %
    % Usage:
    %			symocts = osymset(qA,qB); (calls get_sympairs once per function call)
    %
    %			symocts = osymset(qA,qB,Slist);
    %
    % Dependencies:
    %		get_sympairs.m (optional, required if Slist not supplied)
    %			--allcomb.m
    %
    %		qmult.m
    %
    % Notes:
    %  Could be sped up by doing multiple qA/qB pairs at a time instead of a
    %  single qA/qB pair (i.e. batching/vectorizing approach). Would need to
    %  pay attention to stacking order and perhaps better to output as a cell
    %  instead of an array.
    %--------------------------------------------------------------------------
    %number of symmetry operator pairs
    nsyms = size(Spairs, 1);

    %vertically stack copies of quaternions
    qArep = repmat(qA, nsyms, 1);
    qBrep = repmat(qB, nsyms, 1);

    %unpack pairs
    SAlist = Spairs(:, 1:4);
    SBlist = Spairs(:, 5:8);

    %apply symmetry operators
    qSA = qmult(qArep, SAlist, epsijk);
    qSB = qmult(qBrep, SBlist, epsijk);

    if grainexchangeQ && doublecoverQ
        %apply grain exchange & double cover
        symocts = [ ...
                qSA qSB
                qSA -qSB
                -qSA qSB
                -qSA -qSB
                qSB qSA
                qSB -qSA
                -qSB qSA
                -qSB	-qSA];

    elseif grainexchangeQ && ~doublecoverQ
        symocts = [ ...
                    qSA qSB
                qSB qSA];

    elseif ~grainexchangeQ && doublecoverQ
        symocts = [ ...
                    qSA qSB
                -qSA qSB
                qSA -qSB
                -qSA -qSB];

    elseif ~(grainexchangeQ || doublecoverQ)
        symocts = [ ...
                    qSA qSB];
    end

    %reduce to unique set of octonions
    if uniqueQ
        symocts = uniquetol(round(symocts, 12), 'ByRows', true);
    end

end %osymset

function osets = osymsets(oct, pgnum, usv, grainexchangeQ, doublecoverQ, uniqueQ, epsijk)

    arguments
        oct(:, 8) double {mustBeFinite, mustBeReal, mustBeSqrt2Norm}
        pgnum(1, 1) double {mustBeInteger} = 32 %default to Oh cubic point group
        usv = []
        grainexchangeQ(1, 1) logical {mustBeLogical} = false
        doublecoverQ(1, 1) logical {mustBeLogical} = false
        uniqueQ(1, 1) logical {mustBeLogical} = false
        epsijk(1, 1) double {mustBeInteger} = 1
    end

    % OSYMSETS  Get symmetrically equivalent octonions (osymsets) for each octonion in a list of octonions
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-15
    %
    % Inputs:
    %		data	===	rows of octonions
    %
    %		pgnum ===	point group number (e.g. 32 == cubic)
    %
    %		usv	===	optional, use if "data" is 7D and needs to be projected
    %						back to 8D
    %
    % Outputs:
    %		olist ===	1*size(data,1) cell array containing rows of
    %						unique symmetrically equivalent octonions
    %
    % Dependencies:
    %		Pgnames.mat, PGsymops.mat
    %
    %		osymset.m
    %			--qmult.m
    %
    % Notes:
    %		Adapted a portion from Grain Boundary Octonion function GBdist.m from
    %		Elizabeth Holm's CMU group github page.
    %--------------------------------------------------------------------------

    %% load symmetry operator combinations
    Spairs = get_sympairs(pgnum);

    %% reformat data from 7D Cartesian to 8D Cartesian (if applicable)
    ndatapts = size(oct, 1);

    if size(oct, 2) == 7 && ~isempty(usv)
        oct = proj_up(oct, usv);
    elseif size(oct, 2) == 7
        oct = [oct zeros(size(oct, 1), 1)];
    end

    %% get symmetrically equivalent octonions
    %initialize
    osets = cell(1, ndatapts);

    %unpack quaternions
    qAlist = oct(:, 1:4);
    qBlist = oct(:, 5:8);

    %normalize quaternions
    qAlist = normr(qAlist);
    qBlist = normr(qBlist);

    %loop through quaternion pairs
    parfor i = 1:ndatapts %parfor compatible
        %unpack quaternions
        qA = qAlist(i, :);
        qB = qBlist(i, :);

        %get symmetrically equivalent octonions
        osets{i} = osymset(qA, qB, Spairs, grainexchangeQ, doublecoverQ, uniqueQ, epsijk);
    end

end %osymsets

function [projpts, usv, zeropt] = proj_down(pts, tol, usv, NV)

    arguments
        pts double {mustBeFinite, mustBeReal}
        tol(1, 1) double {mustBeFinite, mustBeReal} = 1e-5
        usv struct = struct.empty
        NV.nforcedim double {mustBeNonnegative, mustBeInteger} = 1
        NV.force(1, 1) logical {mustBeLogical} = false
        NV.zero(1, 1) logical = true
        NV.finaldim double = []
    end

    % PROJ_DOWN  project down by removing null dimensions (i.e. a rotation and translation) via singular value decomposition
    % (SVD)
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Usage:
    %		[projpts,usv] = proj_down(pts);
    %
    %		[projpts,usv] = proj_down(pts,usv);
    %
    % Date:
    %
    % Inputs:
    %   pts - rows of pts to be projected
    %
    %   tol - tolerance for which to consider a dimension "degenerate" in the
    %   square S matrix from SVD.
    %
    %   usv - struct, containing 'U', 'S', 'V' (SVD output), 'avg' (mean of
    %   input points), and possibly 'zeropt' (zeros(1,d) appropriately rotated
    %   and translated via SVD, to preserve circles/spheres, etc.).
    %
    % Outputs:
    %
    % Dependencies:
    %
    % Notes:
    %	Small numerical errors can accumulate by calling proj_down & proj_up
    %	repeatedly with different usv matrices. Shouldn't be an issue if you
    %	keep using the same usv matrices.
    %
    % References:
    %  https://www.mathworks.com/matlabcentral/answers/352830
    %
    %--------------------------------------------------------------------------
    % unpackage name value pairs
    nforcedim = NV.nforcedim;
    forceQ = NV.force;
    finaldim = NV.finaldim;

    % --if zeropt is requested as output, then set zeroQ = true
    if nargout == 3
        zeroQ = true;
    else
        zeroQ = NV.zero;
    end

    %dimensionality
    d = size(pts, 2);

    if ~isempty(finaldim) && forceQ
        nforcedim = d - finaldim;
    end

    if nforcedim >= d
        error(['nforce should be less than d == ' int2str(size(pts, 2))])
    end

    npts = length(pts);

    if ~isempty(usv)
        %unpackage usv
        V = usv.V;
        S = usv.S;
        avg = usv.avg;

        if zeroQ
            zeropt = usv.zeropt;
            Zero = zeros(1, d);
            %projection
            projptstmp = ([Zero; pts] - avg) / V';
            projpts = projptstmp(2:npts, :);
            assert(ismembertol(projptstmp(1, 1:d - nforcedim), zeropt, 1e-6, 'ByRows', true), ...
                ['Zero [' num2str(projptstmp(1, :)) '] did not map back to zeropt [' num2str(zeropt) ']'])
        else
            projpts = (pts - avg) / V';
        end

        %number of degenerate dimensions
        if forceQ
            ndegdim = nforcedim;

        elseif size(S, 1) == size(S, 2)
            ndegdim = sum(abs(diag(S)) < tol);

        else
            %check for columns of all zeros within tolerance
            ndegdim = sum(all(S < tol));
        end

        ndim = size(projpts, 2);

        check_tol = all(abs(projpts(:, ndim - ndegdim + 1:ndim)) < tol, 'all');

        if check_tol
            %remove last column(s)
            projpts = projpts(:, 1:ndim - ndegdim);
        elseif forceQ
            projpts = projpts(:, 1:ndim - ndegdim);
            warning(['Nonzero last column. E.g. ' num2str(pts([1 2],ndim).') ...
                    '. Forcing projection ' int2str(ndegdim) ' dimensions.'])
        elseif ~forceQ
            projpts = pts;
            usv = struct.empty;
            %not sure if I should have a constant, non-zero last column be OK
            if size(pts, 1) > 3
                n = 3;
            else
                n = 1;
            end

            warning(['Nonzero last column. E.g. ' num2str(pts(1:n,ndim).') '. Setting projpts == pts'])
        end

    elseif isempty(usv)
        % make a non-empty struct with no fields
        usv(1) = struct();

        %take average of points
        avg = mean(pts);

        %project to d-1 dimensional space
        if zeroQ
            Zero = zeros(1, size(pts, 2));
            ptstmp = [Zero; pts];
        else
            ptstmp = pts;
        end

        [U, S, V] = svd(ptstmp - avg, 0);

        usv.U = U;
        usv.S = S;
        usv.V = V;
        usv.avg = avg;

        %number of degenerate dimensions
        if forceQ
            ndegdim = nforcedim;

        elseif size(S, 1) == size(S, 2)
            ndegdim = sum(abs(diag(S)) < tol);

        else
            %check for columns of all zeros within tolerance
            ndegdim = sum(all(S < tol));
        end

        if (ndegdim > 0) || forceQ
            %project to lower dimension (i.e. rotation and translation)
            projptstmp = U * S(:, 1:d - ndegdim);

            if zeroQ
                %unpack
                projpts = projptstmp(2:npts, :);
                zeropt = projptstmp(1, :);
                %package
                usv.zeropt = zeropt;
            else
                projpts = projptstmp;
            end

            if (ndegdim == 0) && forceQ
                warning(['ndegdim == 0, tol == ' num2str(tol) ...
                        ', min(diag(S)) == ' num2str(min(diag(S))) ...
                        ', max(diag(S)) == ' num2str(max(diag(S))) ...
                        '. Forcing projection ' int2str(ndegdim) ' dimensions'])
            end

        else
            warning(['ndegdim == 0, tol == ' num2str(tol) ...
                    ', min(diag(S)) == ' num2str(min(diag(S))) ...
                    ', max(diag(S)) == ' num2str(max(diag(S))) ...
                '. Setting projpts == pts'])
            projpts = pts;
            usv = struct.empty;
            zeroQ = false; %override zeroQ
        end

    end

    if zeroQ
        projpts = projpts - zeropt;
    end

end %proj_down

function newpts = proj_up(pts, usv)

    arguments
        pts double
        usv struct {mustContainFields(usv, {'V'})}
    end

    % PROJ_UP  project up (restore null dimensions) using "USV" struct from proj_down.m
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date:
    %
    % Inputs:
    %
    % Outputs:
    %
    % Dependencies:
    %
    %--------------------------------------------------------------------------
    % account for "zeropt" being in usv (see proj_down.m)
    if isfield(usv, 'zeropt')

        if ~isempty(usv.zeropt)
            pts = pts + usv.zeropt;
        end

    end

    V = usv.V;
    avg = usv.avg;

    %lower dimension
    d1 = size(pts, 2);

    %higher dimension
    d2 = size(V, 2);

    ndegdim = d2 - d1; %number of degenerate dimensions

    newpts = padarray(pts, [0 ndegdim], 'post') * V' + avg;

end

function newvertices = projfacet2hyperplane(nvec, vertices)
    % PROJFACET2HYPERPLANE  project facet vertices onto a hyperplane defined by nvec
    % Take a list of vertices defining a facet with vertices on
    % the unit hypersphere and project them onto a tangent hyperplane at a
    % user-defined point, nvec. Useful for computing spherical barycentric
    % coordinates (sphbary.m).
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-03
    %
    % Inputs:
    %	nvec			===	normalized normal to desired tangent hyperplane (i.e.
    %							point on unit hypersphere)
    %
    %	vertices		===	rows of vertices of facet to project from
    %							hypersphere to tangent hyperplane @ nvec
    %
    % Outputs:
    %	newvertices ===	new vertices of facet (projected from hypersphere to
    %							tangent hyperplane @ nvec)
    %
    % Dependencies:
    %	projray2hyperplane.m
    %
    % References:
    %		https://math.stackexchange.com/q/1256236/798661
    %--------------------------------------------------------------------------

    %% project vertices onto hyperplane
    %initialize
    newvertices = zeros(size(vertices));

    nvertices = size(vertices, 1);

    for i = 1:nvertices
        vtx = vertices(i, :); %vertex
        newvertices(i, :) = projray2hyperplane(nvec, vtx);
    end

end %projfacet2hyperplane

function a = projray2hyperplane(nvec, pt)
    % PROJRAY2HYPERPLANE  Project ray (pt) from unit hypersphere to tangent hyperplane at another point (nvec)
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-03
    %
    % Inputs:
    %	nvec	=== normalized normal to hyperplane
    %
    %	pt		=== point to project from hypersphere to hyperplane
    %
    % Intermediates:
    %
    %	d		===	dimensionality
    %
    %	dist	===	distance from origin to hyperplane (same as radius of
    %					hypersphere)
    %
    %	t		===	parametric value that scales pt to the intersection of pt
    %					and tangent hyperplane.
    %
    % Outputs:
    %
    %	a		===	intersection of pt and tangent hyperplane defined by nvec
    %
    % References:
    %		https://math.stackexchange.com/q/1256236/798661
    %
    %		http://comprna.upf.edu/courses/Master_MAT/3_Optimization/U9_Hyperplanes.pdf
    %--------------------------------------------------------------------------

    %dimensionality
    d = length(nvec);

    %distance
    dist = norm(nvec);

    %parametric value
    t = dist / dot(nvec, pt);

    %intersection
    a = -pt * (-t);

end

function [dataProj, facetPts, dataBary, facetIDs, tvals] = ...
        projray2hypersphere(meshpts, facetPtIDs, datanorm, tol, maxnormQ, invmethod)

    arguments
        meshpts double
        facetPtIDs double
        datanorm double
        tol(1, 1) double = 1e-6
        maxnormQ(1, 1) logical = false
        invmethod char {mustBeMember(invmethod, {'mldivide', 'pinv', 'extendedCross'})} = 'mldivide'
    end

    % PROJRAY2HYPERSPHERE  project ray to hypersphere, compute barycentric coordinates, compute intersecting facet
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-06-24
    %
    % Inputs:
    %
    % Outputs:
    %		facetIDs		===	list of facet IDs (i.e. rows of K triangulation)
    %								for which an intersection was found.
    %
    % Dependencies:
    %		numStabBary.m (if using invMethod == 'extendedCross')
    %
    % References:
    %		https://math.stackexchange.com/q/1256236/798661
    %       https://math.stackexchange.com/q/3731002
    %       https://github.com/tomilov/quickhull/issues/8
    %
    % Notes:
    %		To relax the requirement that pts need to be close to on the unit
    %		sphere, then there's a set of lines in projray2hypersphere.m that can
    %		be changed or removed.
    %
    % 		if (t(j) <= 0.1) || (t(j) >= 2.1)
    % 			posQ(j) = 0;
    % 			continue
    % 		end
    %--------------------------------------------------------------------------

    % invmethod = 'mldivide'; %'pinv', 'mldivide', 'extendedCross'

    adjSize = size(facetPtIDs);
    %%
    %setup vectors and matrices for ray projection
    p = cell(adjSize);
    nmatTemp = cell(adjSize);
    nmat = p;
    nvec = p;
    ddet = zeros(adjSize(1), 1);
    dmat = cell(adjSize(1), 1);

    for j = 1:adjSize(1) %loop through facets

        for k = 1:adjSize(2) %loop through vertices of facet
            %package vertices
            p{j, k} = meshpts(facetPtIDs(j, k), :);
        end

        %package vertex matrix
        nmatTemp{j} = vertcat(p{j, :});
        dmat{j} = nmatTemp{j};

        for k = 1:adjSize(2)
            nmat{j, k} = nmatTemp{j};
            nmat{j, k}(:, k) = 1;
        end

        for k = 1:adjSize(2) %loop through dimensions
            nvec{j}(k) = det(nmat{j, k});
        end

        ddet(j) = det(dmat{j});
    end

    %%
    %project ray onto each facet
    a = cell([adjSize(1), 1]);
    lambda = a; %initialize
    posQ = zeros([adjSize(1), 1]);

    warnID1 = 'MATLAB:nearlySingularMatrix';
    warnID2 = 'MATLAB:singularMatrix';
    warnID3 = 'symbolic:mldivide:InconsistentSystem';

    warning('off', warnID1)
    warning('off', warnID2)
    warning('off', warnID3)

    for j = 1:adjSize(1)

        if ddet(j) ~= 0
            t(j) = ddet(j) / dot(datanorm, nvec{j});

            tol2 = 1e-12;

            if (t(j) <= tol2) || (t(j) >= 1 / tol2)
                posQ(j) = 0;
                continue
            end

            a{j} = -datanorm * -t(j);

            switch invmethod
                case 'mldivide'
                    lambda{j} = (dmat{j}' \ a{j}')'; %barycentric coordinates
                    posQ(j) = all(lambda{j} >= -tol) && ~any(lambda{j} == Inf) && (sum(lambda{j}) >= 1 - tol); % && any(c{j} > 0)
                case 'extendedCross'
                    %compute numerically stable barycentric coordinates
                    lambda{j} = numStabBary(nmatTemp{j}, a{j});
                    % 				disp(lambda{j})
                    posQ(j) = all(lambda{j} >= -tol) && (sum(lambda{j}) >= 1 - tol);
            end

        end

    end

    warning('on', warnID1)
    warning('on', warnID2)
    warning('on', warnID3)

    posQtot = sum(posQ);

    %%
    % find facet that intersects ray and compile relevant data for datapoint
    if posQtot > 0

        idList = find(posQ > 0);

        if (length(idList) >= 2) && maxnormQ
            %disp('taking datapoint with largest norm')
            [~, maxpos] = max(vecnorm(vertcat(a{idList})'));
            id = idList(maxpos);
        else
            id = idList;
        end

        dataProj = a{id}; %datapoint projected onto intersecting facet
        facetPts = vertcat(p{id, :}); %each row is a vertex
        dataBary = lambda{id}; %barycentric datapoints (non-generalized)
        facetIDs = id;
        tvals = t(logical(posQ));
        %disp(dataBary{i});

    else
        dataProj = [];
        facetPts = [];
        dataBary = [];
        facetIDs = [];
        tvals = [];
    end

end %projray2hypersphere

function o = qmA2oct(pA, pB, mA, epsijk)

    arguments
        pA(:, 4) double {mustBeFinite, mustBeReal}
        pB(:, 4) double {mustBeFinite, mustBeReal}
        mA(:, 3) double {mustBeFinite, mustBeReal}
        epsijk(1, 1) double = 1
    end

    % QMA2OCT Convert lab coordinates (pA,pB,mA) to octonions (o).
    % Sample frame quaternions of grain A and grain B and sample frame boundary
    % plane normal pointing outward from grain A towards grain B to octonion
    % with BP normal = [0 0 1];
    %--------------------------------------------------------------------------
    % Author(s): Sterling Baird
    % Date: 2020-12-05
    %
    % Inputs:
    %  pA, pB - Quaternions grains A and B in sample reference frame,
    %    resp.
    %  mA - boundary plane normal (sample reference frame) pointing from grain
    %  A to grain B
    %
    % Outputs:
    %   o - octonion, with BP normal = [0 0 1]
    %
    % Usage:
    %  o = qmA2oct(pA,pB,mA)
    %
    % Dependencies:
    %  vecpair2rmat.m
    %  om2qu.m
    %  qu2om.m
    %  qmult.m
    %  qinv.m
    %
    % References:
    %  (1) Francis, T.; Chesser, I.; Singh, S.; Holm, E. A.; De Graef, M. A
    %  Geodesic Octonion Metric for Grain Boundaries. Acta Materialia 2019,
    %  166, 135–147. https://doi.org/10.1016/j.actamat.2018.12.034.

    %--------------------------------------------------------------------------
    npts = size(pB, 1);
    qR = zeros(npts, 4);

    parfor i = 1:npts
        mAtmp = mA(i, :);
        %rotation matrix to go from mAtmp to [0 0 1]
        R = vecpair2rmat(mAtmp, [0 0 1], 1); %not sure why this has to stay in active interpretation for epsijk==1 *and* epsijk==-1
        %convert to quaternion
        qR(i, :) = om2qu(R, epsijk);
    end

    %apply rotation to pA and pB
    qA1 = qmult(qR, pA, epsijk);
    qB1 = qmult(qR, pB, epsijk);

    o = [qA1 qB1];

end

function [c, intfacetIDs] = sphbary_setup(pts, facetPtIDs, data, tol)
    % SPHBARY_SETUP  sphbary coords and intersections for points
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date:
    %
    % Description:
    %
    % Inputs:
    %
    % Outputs:
    %
    % Dependencies:
    %
    %--------------------------------------------------------------------------

    %initialize
    nfacets = size(facetPtIDs, 1);
    posQ = false(nfacets, 1);
    c = zeros(nfacets, size(data, 2));

    %disable warnings
    warnID1 = 'MATLAB:nearlySingularMatrix';
    warnID2 = 'MATLAB:singularMatrix';
    warnID3 = 'symbolic:mldivide:InconsistentSystem';

    warning('off', warnID1)
    warning('off', warnID2)
    warning('off', warnID3)

    %loop through facets
    for i = 1:nfacets
        idtemp = facetPtIDs(i, :);
        vertices = pts(idtemp, :);
        [c(i, :), Pnew] = sphbary(data, vertices);

        %check if spherical barycentric coordinate reqs are met
        posbaryQ = (sum(c(i, :)) >= 1 - tol) && all(c(i, :) >= -tol);

        if posbaryQ

            %enable warnings
            warning('on', warnID1)
            warning('on', warnID2)
            warning('on', warnID3)

            sphbary(data, vertices, Pnew);
            lstwarn = warning('query', 'last');

            if ~strcmp(lstwarn, {warnID1, warnID2, warnID3})
                posQ(i) = true; %position Q, true == intersecting facet
            end

        end

    end

    %enable warnings
    warning('on', warnID1)
    warning('on', warnID2)
    warning('on', warnID3)

    %find intersecting facet IDs
    if sum(posQ) > 0
        intfacetIDs = find(posQ);
        c = c(posQ, :);
    else
        intfacetIDs = [];
        c = [];
    end

end

function [c, Pnew] = sphbary(data, vertices, varargin)
    % SPHBARY  compute spherical barycentric coordinates of a facet or simplex [1]
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-03
    %
    % Inputs:
    %		data		===	datapoint of interest on hypersphere which defines the
    %							tangent hyperplane to the hypersphere
    %
    %		vertices ===	rows of vertices of facet that have been projected to
    %							tangent hyperplane defined by 'a'
    %
    %		Pnew		===	varargin{1}, vertices projected onto tangent
    %							hyperplane
    %
    % Outputs:
    %		c			=== spherical barycentric coordinates of 'a'
    %
    % Dependencies:
    %		projfacet2hyperplane.m
    %			-projray2hyperplane.m
    %
    % References:
    %		[1] T. Langer, A. Belyaev, H.-P. Seidel, Spherical barycentric
    %		coordinates, Proc. Fourth Eurographics Symp. Geom. Process. (2006)
    %		81�88. http://portal.acm.org/citation.cfm?id=1281957.1281968.
    %
    %		https://math.stackexchange.com/q/1256236/798661
    %
    %--------------------------------------------------------------------------
    if nargin == 3
        Pnew = varargin{1};
    end

    if exist('Pnew', 'var') == 0
        Pnew = projfacet2hyperplane(data, vertices);
    end

    c = (Pnew.' \ data.').';

end

function K = sphconvhulln(pts, subhemiQ, dimtype, tol)

    arguments
        pts double {mustBeFinite, mustBeReal}
        subhemiQ logical = true
        dimtype char {mustBeMember(dimtype, {'low', 'high'})} = 'low'
        tol(1, 1) double = 1e-4
    end

    % Compute the "convex hull" on a spherical surface
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-06
    %
    % Example usage: K = sphconvhulln(pts);
    %
    % Inputs:
    %	pts === rows of cospherical points that define a hypersphere surface.
    %
    % Outputs:
    %	K === triangulation of pts without "undercut" facets.
    %
    % Dependencies:
    %  projectfacet2hyperplane.m
    %   -projectray2hyperplane.m
    %
    % Notes:
    %
    %--------------------------------------------------------------------------
    %% setup
    d = size(pts, 2); %dimension

    % check complexity
    if d >= 7 && size(pts, 1) > 400
        warning(['d = ' int2str(d) ' and npts = ' ...
                int2str(size(pts, 1)) '. Initial convex hull input may be too big.'])
        slurmQ = true;

        if ~slurmQ
            m = input('Continue? y/n:', 's');

            if ~strcmp(m, 'y') && ~strcmp(m, 'Y')
                return
            end

        end

    end

    if ~subhemiQ
        K = convhulln(pts);

    elseif subhemiQ

        switch dimtype
            case 'low'
                % probably can handle more points b.c. dimension is lower not sure
                % to what extent it will cause distortions in the triangulation,
                % but will be lessened by not normalizing mean(pts)
                projpts = projfacet2hyperplane(mean(pts), pts);
                a = proj_down(projpts, tol, 'zero', false);
                disp('--delaunayn')
                K = delaunayn(a);

            case 'high'
			%{
				compute convex hull with extra point (faster, 2020-07-31, but
				not sure which can handle more points). Negative scale factor
				removes "undercut" facets. Small positive scale factor keeps
				them (e.g. 0.1)
			%}
                scl = -0.1;
                extrapt = scl * normr(mean(pts)); %assumes points fall on less than a hemisphere

                K = convhulln([pts; extrapt]);

                %remove everything connected to extra point
                npts = size(pts, 1);
                [row, ~] = find(K == npts + 1);
                K(row, :) = [];
        end

    end

end

function pts = sqrt2norm(pts, type, nv)

    arguments
        pts double {mustBeNumeric, mustBeFinite}
        type char {mustBeMember(type, {'oct', 'quat'})} = 'oct'
        nv.warnQ(1, 1) logical = true
    end

    % sqrt2norm  take a set of octonions and give each row norm == sqrt(2) if (norm == 1) || (norm == sqrt(2))
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-27
    %
    % Description: a
    %
    % Inputs:
    %		a -	a
    %
    % Outputs:
    %		b -	b
    %
    % Usage:
    %		a = b(a);
    %
    % Dependencies:
    %		*
    %
    % Notes:
    %		*
    %--------------------------------------------------------------------------
    tol = 1e-2;
    warnQ = nv.warnQ;

    switch type
        case 'quat'
            % if norm(o) == 1 within tolerance, multiply by sqrt(2)
            if all(abs(vecnorm(pts(:, 1:4), 2, 2) - 1 / sqrt(2)) < tol) ...
                && ...
                    all(abs(vecnorm(pts(:, 5:8), 2, 2) - 1 / sqrt(2)) < tol)

                %fine, no warning

            elseif all(abs(vecnorm(pts(:, 1:4), 2, 2) - 1) >= tol) ...
                || ...
                    all(abs(vecnorm(pts(:, 5:8), 2, 2) - 1) >= tol)

                if warnQ
                    warning(['norm(qA) == ' num2str(norm(pts(1, 1:4))) ', ' ...
                            'norm(qB) == ' num2str(norm(pts(1, 5:8))) ...
                        '. norm of 1+ quaternions ~= 0.7071 || 1'])
                end

            end

            pts(:, 1:4) = normr(pts(:, 1:4));
            pts(:, 5:8) = normr(pts(:, 5:8));

        case 'oct'
            ptnm = vecnorm(pts, 2, 2);

            if ~isempty(find(ptnm == 0, 1))
                %            warning('identity octonion(s) present')
                ptnm(ptnm == 0) = [];
            end

            if all(abs(ptnm - 1) < tol)
                %fine, no warning
            elseif any(abs(ptnm - sqrt(2)) >= tol)

                if warnQ
                    warning('norm of octonions ~= 1 || sqrt(2)')
                end

            end

            pts = normr(pts) * sqrt(2);
    end

end

function R = vecpair2rmat(v1, v2, epsijk)

    arguments
        v1(1, 3) double
        v2(1, 3) double
        epsijk(1, 1) double = 1
    end

    % VECPAIR2RMAT  Compute a (non-unique) rotation matrix to go from v1 to v2.
    %  If v1 == v2 or v1 == -v2 within the given numerical precision, then the
    %  identity matrix or -1*identity matrix is given, respectively. Active
    %  rotation matrix, right-handed coordinate system, positive angle of
    %  rotation is counter-clockwise looking from endpoint towards origin,
    %  left-multipled matrix (i.e. R*v1).
    %--------------------------------------------------------------------------
    % Author: Sterling Baird
    %
    % Date: 2020-07-01
    %
    % Inputs:
    %  v1, v2 - two vectors to find the rotation matrix to go from v1->v2
    %
    % Outputs:
    %  R - Rotation matrix such that R*v1 == v2, and R\v2 == v1
    %
    % References
    %	https://math.stackexchange.com/a/897677
    %	https://math.stackexchange.com/a/476311/76513
    %   https://en.wikipedia.org/wiki/Active_and_passive_transformation
    %   https://en.wikipedia.org/wiki/Rotation_matrix
    %
    %   (1) Rowenhorst, D.; Rollett, A. D.; Rohrer, G. S.; Groeber, M.;
    %   Jackson, M.; Konijnenberg, P. J.; De Graef, M. Consistent
    %   Representations of and Conversions between 3D Rotations. Modelling
    %   Simul. Mater. Sci. Eng. 2015, 23 (8), 083501.
    %   https://doi.org/10.1088/0965-0393/23/8/083501.
    %
    % Consider vectorizing
    %--------------------------------------------------------------------------

    precision = 12;
    r = @(a, b) round(a - b, precision);

    isEqual = all(r(v1, v2) == 0); %pointing in same direction
    isOpposite = all(r(v1, -v2) == 0); %pointing in opposite direction
    isParallel = isEqual || isOpposite; %either of previous two cases

    if ~isParallel
        ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        R = eye(3) + ssc(cross(v1, v2)) + ssc(cross(v1, v2))^2 * (1 - dot(v1, v2)) / (norm(cross(v1, v2))^2);

    elseif isEqual
        % same vector
        R = eye(3);

    elseif isOpposite
        % vectors pointing in opposite directions
        R = -eye(3);
    end

    if epsijk == -1 %-1 seems to produce consistent results with section 4.2 of (1)
        R = R.'; %change from active to passive rotation matrix
    end

end

function zm = zeta_min2(o1, o2, epsijk)

    arguments
        o1(:, 8) double {mustBeFinite, mustBeReal}
        o2(:, 8) double {mustBeFinite, mustBeReal}
        epsijk(1, 1) double = 1
    end

    % ZETA_MIN  Alternative version of CMU group function zeta_min(), vectorized by Sterling Baird
    %--------------------------------------------------------------------------
    % Date: 2020-07-27
    %
    % Inputs:
    %		(o1,o2)	-	lists of octonions
    %
    % Outputs:
    %		zm			-	list of minimized zeta angles
    %
    % Usage:
    %		zm = zeta_min2(o1,o2);
    %
    % Dependencies:
    %		*
    %
    % Notes:
    %		*
    %--------------------------------------------------------------------------

    %quaternion dot products = cos(omega/2), omega is misorientation angle

    %unpack quaternions
    qA = o1(:, 1:4);
    qB = o1(:, 5:8);
    qC = o2(:, 1:4);
    qD = o2(:, 5:8);

    %dot() applies conj, so be careful that inputs are real. Alternative: sum(qA.*qC);
    qdot_AC = dot(qA, qC, 2);
    qdot_BD = dot(qB, qD, 2);

    mu_num1 = qA(:, 4) .* qC(:, 1) - qC(:, 4) .* qA(:, 1) + qB(:, 4) .* qD(:, 1) - qD(:, 4) .* qB(:, 1);
    crossAC = cross(qA(:, 2:4), qC(:, 2:4), 2);
    crossBD = cross(qB(:, 2:4), qD(:, 2:4), 2);

    mu_arg = (mu_num1 + epsijk * crossAC(:, 3) + epsijk * crossBD(:, 3)) ./ (qdot_AC + qdot_BD);

    mu_arg(isnan(mu_arg)) = 0;

    mu = 2 * atan(mu_arg);

    % tol = 1e-6;
    % mu(abs(mu - pi) < tol) = pi;
    % mu(abs(mu + pi) < tol) = -pi;

    if any(~isfinite(mu))
        1 + 1;
    end

    % shift negative values
    zm = mu;
    negIDs = mu < 0;
    zm(negIDs) = zm(negIDs) + 2 * pi;

end

function q = ax2qu(ax, epsijk)

    arguments
        ax(:, 4)
        epsijk = 1
    end

    % AX2QU  from axis-angle pair to quaternions, first 3 elements:axis, last element: angle
    %vectorized by SGB 2020-08-15
    npts = size(ax, 1);
    q = zeros(npts, 4);

    thr = 1e-10;
    ids = abs(ax(:, 4)) < thr;

    if any(ids)
        nids = sum(ids);
        q(ids, :) = repmat([1.0, 0.0, 0.0, 0.0], nids, 1);
    end

    if any(~ids)
        c = cos(ax(~ids, 4) * 0.5);
        s = sin(ax(~ids, 4) * 0.5);
        q(~ids, :) = [c, ax(~ids, 1) .* s, ax(~ids, 2) .* s, ax(~ids, 3) .* s];
        q(~ids, 2:4) = -epsijk * q(~ids, 2:4);
    end

    % set values very close to 0 as 0

    q(abs(q) < thr) = 0;
end

% from cubochoric to homochoric

function q = cu2ho(xyz, printQ)

    arguments
        xyz
        printQ(1, 1) logical = false
    end

    % rs = sqrt(sum(xyz.*xyz));
    R1 = (3 * pi / 4)^(1/3);
    aa = ((pi^5) / 6)^(1/6);
    bb = aa / 2;
    sc = aa / (pi^(2/3));
    prek = R1 * (2^0.25) / bb;

    % ierr = 0;
    if (max(abs(xyz)) > ((pi^(2/3) / 2.0) + 1e-8))
        q = [0.0, 0.0, 0.0];
        %   ierr = 1;
        return
    end

    % determine which pyramid pair the point lies in and copy coordinates in correct order (see paper)
    if printQ
        disp('GetPyramid')
    end

    p = GetPyramid(xyz);

    if (p == 1) || (p == 2)
        sXYZ = xyz;
    elseif (p == 3) || (p == 4)
        sXYZ = [xyz(2), xyz(3), xyz(1)];
    else
        sXYZ = [xyz(3), xyz(1), xyz(2)];
    end

    % scale by grid parameter ratio sc
    XYZ = sc * sXYZ;

    % transform to the sphere grid via the curved square, and intercept the zero point
    if (max(abs(XYZ)) == 0.0)
        LamXYZ = [0.0, 0.0, 0.0];
    else
        % intercept all the points along the z-axis
        if (max(abs(XYZ(1:2))) == 0.0)
            LamXYZ = [0.0, 0.0, sqrt(6 / pi) * XYZ(3)];
        else % this is a general grid point

            if (abs(XYZ(2)) <= abs(XYZ(1)))
                c = cos((pi / 12) * XYZ(2) / XYZ(1));
                s = sin((pi / 12) * XYZ(2) / XYZ(1));
                temp = prek * XYZ(1) / sqrt(sqrt(2) - c);
                T1 = (sqrt(2) * c - 1.0) * temp;
                T2 = sqrt(2) * s * temp;
            else
                c = cos((pi / 12) * XYZ(1) / XYZ(2));
                s = sin((pi / 12) * XYZ(1) / XYZ(2));
                temp = prek * XYZ(2) / sqrt(sqrt(2) - c);
                T1 = sqrt(2) * s * temp;
                T2 = (sqrt(2) * c - 1.0) * temp;
            end

            % transform to sphere grid (inverse Lambert)
            % [note that there is no need to worry about dividing by zero, since XYZ(3) can not become zero]
            c = T1^2 + T2^2;
            s = pi * c / (24.0 * XYZ(3)^2);
            c = sqrt(pi) * c / sqrt(24) / XYZ(3);
            q = sqrt(1.0 - s);
            LamXYZ = [T1 * q, T2 * q, sqrt(6 / pi) * XYZ(3) - c];
        end

    end

    % reverse the coordinates back to the regular order according to the original pyramid number
    if (p == 1) || (p == 2)
        q = LamXYZ;
    elseif (p == 3) || (p == 4)
        q = [LamXYZ(3), LamXYZ(1), LamXYZ(2)];
    else
        q = [LamXYZ(2), LamXYZ(3), LamXYZ(1)];
    end

    % set values very close to 0 as 0
    thr = 1e-7;

    if (abs(q(1)) - 0) < thr
        q(1) = 0;
    elseif (abs(q(2)) - 0) < thr
        q(2) = 0;
    elseif (abs(q(3)) - 0) < thr
        q(3) = 0;
    end

end

% from cubochoric to quaternions

function q = cu2qu(c, printQ)

    arguments
        c
        printQ(1, 1) logical = false
    end

    if printQ
        disp('cu2ho')
    end

    h = cu2ho(c, printQ);

    if printQ
        disp('ho2qu')
    end

    q = ho2qu(h, printQ);

    % set values very close to 0 as 0
    thr = 1e-10;

    if (abs(q(1)) - 0) < thr
        q(1) = 0;
    elseif (abs(q(2)) - 0) < thr
        q(2) = 0;
    elseif (abs(q(3)) - 0) < thr
        q(3) = 0;
    elseif (abs(q(4)) - 0) < thr
        q(4) = 0;
    end

end

function q = GetPyramid(xyz)

    next = 1;

    if ((abs(xyz(1)) <= xyz(3)) && (abs(xyz(2)) <= xyz(3)))
        p = 1; % pyramid 1
        next = 0;
    end

    if (next)

        if ((abs(xyz(1)) <= -xyz(3)) && (abs(xyz(2)) <= -xyz(3)))
            p = 2; % pyramid 2
            next = 0;
        end

    end

    if (next)

        if ((abs(xyz(3)) <= xyz(1)) && (abs(xyz(2)) <= xyz(1)))
            p = 3; % pyramid 3
            next = 0;
        end

    end

    if (next)

        if ((abs(xyz(3)) <= -xyz(1)) && (abs(xyz(2)) <= -xyz(1)))
            p = 4; % pyramid 4
            next = 0;
        end

    end

    if (next)

        if ((abs(xyz(1)) <= xyz(2)) && (abs(xyz(3)) <= xyz(2)))
            p = 5; % pyramid 5
            next = 0;
        end

    end

    if (next)

        if ((abs(xyz(1)) <= -xyz(2)) && (abs(xyz(3)) <= -xyz(2)))
            p = 6; % pyramid 6
            %   next = 0;
        end

    end

    q = p;

end

% from homochoric to axis-angle pair

function q = ho2ax(h)

    thr = 1e-8;
    tfit = [+1.0000000000018852, -0.5000000002194847, ...
            -0.024999992127593126, - 0.003928701544781374, ...
            -0.0008152701535450438, - 0.0002009500426119712, ...
            -0.00002397986776071756, - 0.00008202868926605841, ...
            +0.00012448715042090092, - 0.0001749114214822577, ...
            +0.0001703481934140054, - 0.00012062065004116828, ...
            +0.000059719705868660826, - 0.00001980756723965647, ...
            +0.000003953714684212874, - 0.00000036555001439719544];

    % normalize xyz and store the magnitude
    hmag = sum(h .* h);

    if (hmag == 0.0)
        q = [0.0, 0.0, 1.0, 0.0];
    else
        hm = hmag;
        hn = h / sqrt(hmag);

        %convert the magnitude to the rotation angle
        s = tfit(1) + tfit(2) * hmag;

        for i = 3:16
            hm = hm * hmag;
            s = s + tfit(i) * hm;
        end

        s = 2.0 * acos(s);

        if (abs(s - pi) < thr)
            q = [hn(1), hn(2), hn(3), pi];
        else
            q = [hn(1), hn(2), hn(3), s];
        end

    end

    % set values very close to 0 as 0
    if (abs(q(1)) - 0) < thr
        q(1) = 0;
    elseif (abs(q(2)) - 0) < thr
        q(2) = 0;
    elseif (abs(q(3)) - 0) < thr
        q(3) = 0;
    end

end

% from homochoric to quaternions

function q = ho2qu(h, printQ)

    arguments
        h
        printQ(1, 1) logical = false
    end

    if printQ
        disp('ho2ax')
    end

    ax = ho2ax(h);

    if printQ
        disp('ax2qu')
    end

    q = ax2qu(ax);

    % set values very close to 0 as 0
    thr = 1e-8;

    if (abs(q(1)) - 0) < thr
        q(1) = 0;
    elseif (abs(q(2)) - 0) < thr
        q(2) = 0;
    elseif (abs(q(3)) - 0) < thr
        q(3) = 0;
    elseif (abs(q(4)) - 0) < thr
        q(4) = 0;
    end

end

function q = om2ax(om, epsijk)

    arguments
        om(3, 3) double
        epsijk(1, 1) double = 1
    end

    % OM2AX  from rotation matrix to axis-angle pair

    % global epsijk

    q = zeros(1, 4);
    thr = 1e-10;

    % first get the rotation angle
    t = 0.5 * (om(1, 1) + om(2, 2) + om(3, 3) - 1.0);

    if t > 1.0
        t = 1.0;
    elseif t <- 1.0
        t = -1.0;
    end

    q(4) = acos(t);

    if (q(4) == 0.0)
        q(1:3) = [0.0, 0.0, 1.0];
        return
    else
        % find eigenvectors V and eigenvalues W of the rotation matrix
        [V, W] = eig(om);

        for i = 1:3

            if ((abs(real(W(i, i)) - 1.0) < thr) && (abs(imag(W(i, i))) < thr))
                q(1:3) = V(:, i);
                % change sign
                % matlab chose descending order for same eigenvector value
                % eg. if v1=-v2, it would choose v1=-1 and v2=1,
                % need to change to v1=1 and v2=-1.
                if q(1) < 0
                    q(1:3) = -q(1:3);
                end

                if ((om(2, 3) - om(3, 2)) ~= 0.0)
                    s = sign(-epsijk * (om(2, 3) - om(3, 2)));
                    q(1) = s * abs(q(1));
                end

                if ((om(3, 1) - om(1, 3)) ~= 0.0)
                    s = sign(-epsijk * (om(3, 1) - om(1, 3)));
                    q(2) = s * abs(q(2));
                end

                if ((om(1, 2) - om(2, 1)) ~= 0.0)
                    s = sign(-epsijk * (om(1, 2) - om(2, 1)));
                    q(3) = s * abs(q(3));
                end

                return
            end

        end

    end

end

function q = om2qu(omIn, epsijk)

    arguments
        omIn
        epsijk(1, 1) double = 1
    end

    % OM2QU  from rotation matrix to quaternions

    % global epsijk

    if isnumeric(omIn)
        npts = 1;
    else
        npts = length(omIn);
    end

    qout = zeros(npts, 4);

    for i = 1:npts

        if isnumeric(omIn)
            om = omIn;
        elseif iscell(omIn)
            om = omIn{i};
        end

        thr = 1e-10;

        s = 1 + om(1, 1) + om(2, 2) + om(3, 3);

        if abs(s) < thr
            s = 0;
        end

        q0 = 0.5 * sqrt(s);

        s = 1 + om(1, 1) - om(2, 2) - om(3, 3);

        if abs(s) < thr
            s = 0;
        end

        q1 = 0.5 * sqrt(s);

        s = 1 - om(1, 1) + om(2, 2) - om(3, 3);

        if abs(s) < thr
            s = 0;
        end

        q2 = 0.5 * sqrt(s);

        s = 1 - om(1, 1) - om(2, 2) + om(3, 3);

        if abs(s) < thr
            s = 0;
        end

        q3 = 0.5 * sqrt(s);

        % verify the signs (q0 always positive)
        if epsijk == 1

            if om(3, 2) < om(2, 3)
                q1 = -epsijk * q1;
            end

            if om(1, 3) < om(3, 1)
                q2 = -epsijk * q2;
            end

            if om(2, 1) < om(1, 2)
                q3 = -epsijk * q3;
            end

        end

        % normalization
        magq = sqrt(q0^2 + q1^2 + q2^2 + q3^2);

        if magq ~= 0.0
            q = [q0, q1, q2, q3] / magq;
        end

        % ! we need to do a quick test here to make sure that the
        % ! sign of the vector part is the same as that of the
        % ! corresponding vector in the axis-angle representation;
        % ! these two can end up being different, presumably due to rounding
        % ! issues, but this needs to be further analyzed...
        % ! This adds a little bit of computation overhead but for now it
        % ! is the easiest way to make sure the signs are correct.

        oax = om2ax(om, epsijk);

        if (oax(1) * q(2) < 0.0)
            q(2) = -q(2);
        end

        if (oax(2) * q(3) < 0.0)
            q(3) = -q(3);
        end

        if (oax(3) * q(4) < 0.0)
            q(4) = -q(4);
        end

        %package
        qout(i, :) = q;
    end

end

function om = qu2om(qq, epsijk, method)

    arguments
        qq(:, 4) double
        epsijk(1, 1) double = 1
        method = 'new'
    end

    % QU2OM  from quaternions to rotation matrix

    npts = size(qq, 1);
    % global epsijk
    switch method
        case 'new'
            qbar = qq(:, 1) .* qq(:, 1) - (qq(:, 2) .* qq(:, 2) + qq(:, 3) .* qq(:, 3) + qq(:, 4) .* qq(:, 4));

            om = zeros(3, 3, npts);

            fn = @(x) reshape(x, 1, 1, npts);
            om(1, 1, :) = fn(qbar + 2.0 * qq(:, 2) .* qq(:, 2));
            om(1, 2, :) = fn(2.0 * (qq(:, 2) .* qq(:, 3) - qq(:, 1) .* qq(:, 4)));
            om(1, 3, :) = fn(2.0 * (qq(:, 2) .* qq(:, 4) + qq(:, 1) .* qq(:, 3)));
            om(2, 1, :) = fn(2.0 * (qq(:, 3) .* qq(:, 2) + qq(:, 1) .* qq(:, 4)));
            om(2, 2, :) = fn(qbar + 2.0 * qq(:, 3) .* qq(:, 3));
            om(2, 3, :) = fn(2.0 * (qq(:, 3) .* qq(:, 4) - qq(:, 1) .* qq(:, 2)));
            om(3, 1, :) = fn(2.0 * (qq(:, 4) .* qq(:, 2) - qq(:, 1) .* qq(:, 3)));
            om(3, 2, :) = fn(2.0 * (qq(:, 4) .* qq(:, 3) + qq(:, 1) .* qq(:, 2)));
            om(3, 3, :) = fn(qbar + 2.0 * qq(:, 4) .* qq(:, 4));

            thr = eps;
            om(abs(om) < thr) = 0;

            if epsijk ~= 1
                om = permute(om, [2 1 3]);
            end

        case 'original'
            qbar = qq(1) * qq(1) - (qq(2) * qq(2) + qq(3) * qq(3) + qq(4) * qq(4));

            om(1, 1) = qbar + 2.0 * qq(2) * qq(2); %om1
            om(2, 2) = qbar + 2.0 * qq(3) * qq(3); %om5
            om(3, 3) = qbar + 2.0 * qq(4) * qq(4); %om9

            om(1, 2) = 2.0 * (qq(2) * qq(3) - qq(1) * qq(4)); %om2
            om(2, 3) = 2.0 * (qq(3) * qq(4) - qq(1) * qq(2)); %om6
            om(3, 1) = 2.0 * (qq(4) * qq(2) - qq(1) * qq(3)); %om7
            om(2, 1) = 2.0 * (qq(3) * qq(2) + qq(1) * qq(4)); %om4
            om(3, 2) = 2.0 * (qq(4) * qq(3) + qq(1) * qq(2)); %om8
            om(1, 3) = 2.0 * (qq(2) * qq(4) + qq(1) * qq(3)); %om3

            if (epsijk ~= 1)
                om = transpose(om);
            end

            thr = 1e-8;

            for i = 1:3

                for j = 1:3

                    if (abs(om(i, j)) < thr)
                        om(i, j) = 0.0;
                    end

                end

            end

    end

end

function out = qinv(q)

    arguments
        q double
    end

    %vectorized by SGB 2020-08-15
    %%% Lp(r) = vec(prp*) active. For passive, implement Lp*(r) = vec(p*rp)
    % also see https://dx.doi.org/10.1088/0965-0393/23/8/083501 eq. 24 in
    % regards to Lp(r)

    q0 = q(:, 1);
    q = q(:, 2:4);
    out = [q0 -q];

end

function out = qmult(pp, qq, epsijk)

    arguments
        pp(:, 4) double {mustBeReal, mustBeFinite}
        qq(:, 4) double {mustBeReal, mustBeFinite}
        epsijk = 1
    end

    % multiply lists of quaternion pairs (input: rows of quaternions),
    % vectorized by SGB 2020-07-27
    p = pp(:, 2:4); q = qq(:, 2:4);
    qr = pp(:, 1) .* qq(:, 1) - dot(p, q, 2);
    qi = pp(:, 1) .* q + qq(:, 1) .* p + epsijk * cross(p, q, 2);

    out = [qr qi];

end

function [qm, nA, y] = get_olmsted_ni_data()
    qm = [
        0 0.948683300167071 0.316227759667166 0
        0.948683297571594 0.316227767453598 0 0
        0 0.577350267951537 0.577350271803421 0.577350267813920
        0 -0.816496581900497 -0.408248283238622 -0.408248295743562
        -0.577350267813920 0.577350271803421 0.577350267951537 0
        0 0.894427182880770 0.447213611738249 0
        0 0.980580674639535 0.196116140395112 0
        0.980580674323291 0.196116141976329 0 0
        0 -0.912870925491948 -0.365148377716096 -0.182574192159729
        -0.894427191645603 -0.447213594208583 0.000000003880661 0
        0.408248286220753 -0.816496578931194 -0.408248298700037 0
        0.912870925491948 -0.365148377716096 -0.182574192159729 0
        0.408248298700037 0.816496578931194 0.408248286220753 0
        -0.774596668324746 0.516397784326625 0.258198887100898 0.258198885479032
        0.632455533156445 0.632455532355850 0.316227765161674 0.316227763982115
        0 0.857492918167444 0.514495768002693 0
        0.857492919604163 0.514495765608161 0 0
        0.288675137521039 0.866025405287300 0.288675130466852 0.288675131287965
        0 0.666666675251644 0.666666660164754 0.333333329167204
        -0.408248285206660 -0.816496583060147 -0.408248291456225 0.000000000418802
        0 -0.942809043055013 -0.235702255875100 -0.235702259024135
        -0.707106779612030 0.471404520799178 0.471404520799178 0.235702265086485
        0.166666674317966 -0.833333334787778 -0.499999994989542 -0.166666666774521
        -0.333333330387435 0.666666672737162 0.666666662069120 0
        0 0.771516756494377 0.617213391055954 0.154303351712157
        0.801783725639041 0.534522487575102 0.267261234706613 0
        -0.577350268858319 0.577350275117051 0.577350263593508 0
        0 -0.801783725639041 -0.534522487575102 -0.267261234706613
        -0.771516756494377 -0.617213391055954 -0.154303351712157 0
        0.577350263593508 -0.577350275117051 -0.577350268858319 0
        0.981980505361639 -0.109108949072803 -0.109108943895178 -0.109108948688952
        0.188982240143620 0.566946710294535 0.566946704518728 0.566946712515257
        0 -0.904534031595092 -0.301511355883906 -0.301511339686217
        0 0.639602149757223 0.639602145352570 0.426401437247025
        0.904534031595092 -0.301511339686217 -0.301511355883906 0
        0 0.989949493510053 0.141421357295107 0
        0.948683298606537 0 0 -0.316227764348769
        -0.387298331188661 0.903696115766523 0.129099445068834 0.129099443414359
        0.989949492443151 0.141421364763421 0 0
        0 0.832050294073944 0.554700196621079 0
        0 0.919145032908829 0.393919291834036 0
        0.919145031813435 0.393919294389954 0 0
        0 0.861640435364518 0.492365966086633 0.123091492737805
        0.904534028759875 -0.301511352534824 -0.301511351540951 0
        0.408248300952997 0.816496576681248 0.408248288467685 0
        0.957427109071476 -0.261116479282734 -0.087038832367706 -0.087038823073916
        -0.288675130232995 -0.866025402928020 -0.288675137325062 -0.288675138795638
        -0.522232965698035 0.696310623706196 0.348155312632244 0.348155314677326
        0.852802866750713 0.426401432798962 0.213200717407906 0.213200709814745
        0 -0.727606878247232 -0.485071240616695 -0.485071254821287
        0 -0.685994337241436 -0.514495755088474 -0.514495760204711
        0 0.970142502121177 0.242535617132954 0
        0.727606878247232 -0.485071254821287 -0.485071240616694 0
        0 0.717137165764438 0.597614307820115 0.358568577217851
        0.894427191756160 -0.000000006208917 -0.447213593987470 0
        0.267261246871213 0.801783723387227 0.534522484870524 0
        -0.597614305636337 0.717137165742673 0.358568580901010 0
        0.801783725014918 0.534522488841321 0.267261234046544 0
        0.358568586723731 -0.717137170076859 -0.597614296941682 0
        0.894427189262067 0.447213598975655 -0.000000002691218 0
        0 0.813733472837210 0.581238191436432 0
        0.813733471766573 0.581238192935324 0 0
        0 0.688247204012338 0.688247198794377 0.229415735120529
        0 -0.973328526489955 -0.162221427587113 -0.162221416442132
        -0.229415729279848 0.688247210011710 0.688247194741899 0
        0 -0.792593921867616 -0.566138517742786 -0.226455412270304
        -0.784464543605915 -0.588348403541864 -0.196116128543533 0
        0.577350262906886 -0.577350270440520 -0.577350274221471 0
        0 0.784464543605915 0.588348403541864 0.196116128543533
        0.792593921867616 0.566138517742786 0.226455412270303 0
        -0.577350274221471 0.577350270440520 0.577350262906886 0
        -0.960768923000843 0.160128149142441 0.160128152481452 0.160128158769448
        0.277350097522609 0.554700194523116 0.554700196806196 0.554700197641378
        0 0.993883735904950 0.110431514992866 0
        0.993883733972952 0.110431532380848 0 0
        0 -0.872871563406380 -0.436435781562180 -0.218217878205962
        -0.534522476668145 -0.801783727140108 -0.267261252017329 0
        0.816496588828049 -0.408248286543283 -0.408248278583796 0
        0.872871563042085 -0.218217898067548 -0.436435772359976 0
        0.408248293874093 0.816496577507801 0.408248293893483 0
        -0.925820102073453 0.308606700460995 0.154303346410528 0.154303338634623
        0.377964467373192 0.755928950808736 0.377964469988790 0.377964472085137
        0.516397778914041 0.774596668443589 0.258198887782670 0.258198895265897
        -0.816496579549551 0.408248290900646 0.408248292783431 -0.000000003853564
        0 0.737864781795359 0.527046281046556 0.421637025676211
        0.843274042330571 0.527046276844005 -0.105409257640543 0
        -0.333333328832231 0.666666663108233 0.666666672475651 0
        0.903696115203099 0.387298332154532 0.129099441023161 0.129099448506387
        -0.408248286906601 0.816496583543596 0.408248288789386 0.000000000140481
        0.833333331554116 -0.500000000653941 -0.166666668445883 -0.166666671821710
        0.527046278945281 0.632455528045131 0.527046278945281 0.210818511390771
        0.948683296080799 -0.000000002101276 -0.316227771925983 -0.000000002101276
        0.166666671821710 0.833333331554116 0.500000000653941 -0.166666668445883
        0 0.843274042330571 0.527046276844005 0.105409257640543
        -0.737864781795359 0.527046281046556 0.421637025676211 0
        0.666666672475651 0.666666663108233 0.333333328832231 0
        0.645497223510346 -0.645497225225460 -0.387298334792253 0.129099444359047
        -0.447213600666561 -0.596284792277743 -0.596284792277743 -0.298142396138871
        0.894427188416614 -0.298142400444374 -0.298142400444374 -0.149071200222187
        -0.948683299369664 -0.210818508039591 -0.210818508039591 -0.105409254019795
        0.316227762059386 -0.632455532913110 -0.632455532913110 -0.316227766456555
        0.670820394664150 0.670820394664150 0.223606793507338 0.223606793507338
        0 0.800000007589466 0.599999989880711 0
        0.424264072655926 0.848528136307565 0.282842707066162 0.141421361919987
        -0.707106778519293 0.565685424178679 0.424264074184766 0
        0.900000002036494 -0.299999996310820 0.100000003806633 -0.299999996310820
        0.282842707066162 0.848528136307565 0.424264072655926 0.141421361919987
        0.948683298815881 0.316227763720737 -0.000000005752621 -0.000000001219727
        0 0.707106778519293 0.565685424178679 0.424264074184766
        0.800000002904363 0.599999996127516 -0.000000011686254 0
        -0.424264057762572 0.565685429827187 0.707106783853803 0
        0 0.874157276519123 0.485642930462977 0
        0.874157276248418 0.485642930950247 0 0
        0 -0.962250448705994 -0.192450081347409 -0.192450097829253
        -0.471404516568252 -0.707106785625854 -0.471404516766169 -0.235702263572879
        0.816496583033894 -0.408248282019401 0.000000001253501 -0.408248294695989
        0 0.680413820500869 0.680413814458744 0.272165526775735
        -0.471404516766169 0.707106785625854 0.471404516568252 0.235702263572880
        0.866025405508607 0.288675130596311 0.288675132267718 0.288675135747905
        0.866025403204854 -0.481125225176152 -0.096225042275991 -0.096225048412825
        -0.500000001003870 -0.833333333768230 -0.166666666052446 -0.166666662094796
        0.952579346757920 -0.136082756636102 0.000000005849733 -0.272165522203995
        0.166666657087664 0.833333335376850 0.499999999275700 0.166666668200985
        0 0.952579344260516 0.272165528969358 0.136082760587208
        0.962250450855370 -0.192450082847319 -0.192450085582461 0
        0.235702253064329 0.942809042411500 0.235702264408955 0
        -0.288675139123406 0.866025403871734 0.288675129902245 0.288675134496900
        0.942809039978450 0.235702262639806 0.235702264565681 0.000000006043862
        -0.816496581067106 0.544331054527627 0.136082760312359 0.136082763524028
        0.577350268992512 0.769800361447263 0.192450087075288 0.192450082864756
        0.952579343187462 -0.000000005001163 -0.136082763439572 -0.272165531298864
        0.962250446758467 -0.192450093021129 -0.192450095893168 0
        0 -0.953462587590819 -0.286038780861493 -0.095346263208841
        0.301511337756518 -0.904534034262755 -0.301511349810614 0
        0.948683299787861 0.316227760804795 -0.000000004644379 0
        -0.741619848107824 0.606779875339475 0.202259965001612 0.202259957315050
        0.670820393915188 0.670820391868726 0.223606802381152 0.223606795266687
        0.953462588937245 -0.286038775607732 -0.095346265505867 0
        0.301511345552842 0.904534036137430 0.301511336390267 0
        0 0.749268645059556 0.655610072251617 0.093658585793540
        0.811107102130947 0.486664267127736 0.324442845465517 0
        -0.577350273386868 0.577350261334706 0.577350272847303 0
        0 -0.811107102130947 -0.486664267127736 -0.324442845465517
        -0.749268645059556 -0.655610072251617 -0.093658585793540 0
        0.577350272847303 -0.577350261334706 -0.577350273386868 0
        0.993399267611945 -0.066226622167489 -0.066226619210995 -0.066226614983838
        0.114707868553343 0.573539337059962 0.573539341551377 0.573539325094263
        0 0.928476690782124 0.371390676611943 0
        0 0.995893206728755 0.090535743170954 0
        0.995893207586201 0.090535733739044 0 0
        0 0.964763821767187 0.263117403850776 0
        0.980580676527078 0 0 -0.196116130957393
        0 0.894427193806550 0.447213589886690 0
        0 0.789352223846899 0.613940605195613 0
        0.964763822493052 0.263117401189269 0 0
        0.789352210704815 0.613940622092578 0 0
        0 0.696310618597245 0.696310628500968 0.174077658145175
        -0.452267017254536 -0.753778357278241 -0.452267023519456 -0.150755671997617
        0.866025403972169 -0.288675134147446 -0.288675135064955 -0.288675134008846
        0 -0.870388277890326 -0.348155319575241 -0.348155308967958
        0 -0.984731927447613 -0.123091491129882 -0.123091493925171
        -0.452267023519456 -0.753778357278241 -0.452267017254536 -0.150755671997617
        0.816496577057503 -0.408248295446840 -0.000000002939850 -0.408248293221331
        0 0.615457452585078 0.615457454537983 0.492365967255163
        -0.174077658145175 0.696310628500969 0.696310618597245 0
        0.870388282104869 -0.348155313608394 -0.348155304398448 0
        0 -0.845154252121480 -0.507092559738786 -0.169030843275855
        0.267261242012945 -0.801783722955709 -0.534522487946935 0
        0.948683299388764 0.316227762002087 0.000000008515254 0
        -0.169030852936886 -0.845154258073266 -0.507092546598800 0
        0.948683295867640 -0.000000000203511 -0.316227772565459 0
        -0.845154252566424 0.507092553704450 0.169030859154148 0
        0.534522487243417 0.801783724561420 0.267261238602846 0
        -0.645497222094107 -0.645497222758091 -0.387298338414243 -0.129099452911111
        0.763762617747681 -0.545544729281514 -0.327326826976617 -0.109108938340493
        -0.552052448554671 -0.759072116121051 -0.345032776084102 0.000000002325056
        0.816496579831817 -0.408248292011731 -0.408248291107812 0.000000013649707
        -0.408248287556758 -0.816496578935353 -0.408248297355713 -0.000000000491209
        0.897085229727995 -0.276026222431027 -0.345032774000279 0.000000003487112
        0.979957885397199 0.178174163495740 -0.089087093961156 0
        -0.105409252397131 0.843274049411578 0.527046266563077 0
        0.986013296110267 -0.140859042980587 -0.084515434457470 -0.028171816463007
        0.166666673014634 0.833333336029861 0.499999992910366 0.166666668104965
        0.890870803965847 0.445435406252094 0.089087089402880 0
        -0.421637028732181 0.737864784879197 0.527046274284407 0
        -0.516397784671447 0.774596664233338 0.258198887909948 0.258198896254561
        0.828078669318008 0.483045897322649 0.276026221891845 -0.069006545547418
        0.872871558663120 0.436435787567972 -0.000000007163111 -0.218217885167415
        -0.408248292980575 0.816496578339495 0.408248293123612 0.000000002173567
        0 0.979957885397199 0.178174163495740 0.089087093961156
        -0.845154259642273 0.507092540914963 0.169030862143363 0
        0.527046266563077 0.843274049411578 0.105409252397131 0
        -0.833333336029861 0.499999992910366 0.166666668104965 0.166666673014634
        0.534522479482329 0.801783727881475 0.267261244164859 0.000000010050198
        0.986013296110267 0.028171816463007 -0.084515434457470 -0.140859042980587
        -0.000000004054657 0.948683296926543 0.316227769388749 0.000000002890705
        0.903696114614905 -0.387298333731615 -0.129099445092482 -0.129099443823172
        -0.414039335411543 -0.759072114768315 -0.483045889969487 -0.138013120741150
        -0.545544727323501 -0.763762617312145 -0.327326830131808 0.109108941713738
        0.816496578435797 -0.408248289472107 -0.408248296439477 0.000000003507143
        0 -0.890870803965847 -0.445435406252094 -0.089087089402880
        -0.845154255307662 -0.507092553773183 -0.169030845241758 0
        0.527046274284407 -0.737864784879197 -0.421637028732181 0
        -0.169030854316980 -0.845154255441746 -0.507092550524635 0
        0.979957886084194 -0.089087082170664 -0.178174165612513 0
        0.986013297718991 -0.084515425867247 0.140859038556939 -0.028171808046566
        -0.000000002668501 -0.801783726887320 -0.534522482048728 -0.267261242014525
        -0.545544725058454 -0.763762615158424 -0.327326836564055 -0.109108948818275
        0.828078671930105 -0.483045891965175 -0.276026221909518 0.069006551633882
        -0.169030848843768 -0.845154256261921 -0.507092550982080 0
        0.890870807996287 0.089087079571544 -0.445435400157481 0
        0.872871561685108 0.436435781348806 -0.000000002639841 0.218217885517795
        0.414039333662615 -0.759072116661881 -0.483045891987445 -0.138013108510471
        0 0.910366473688846 0.413802952603453 0
        0.910366478550112 0.413802941908667 0 0
        0 0.986393924174848 0.164398985249130 0
        0 -0.966987555641265 -0.241746894314021 -0.080582295353615
        -0.507092546270047 -0.845154257382414 -0.169030857377403 0
        0.858116334300252 -0.476731287328517 -0.190692518178328 0
        -0.944911181374135 -0.188982238797995 -0.188982240135976 0.188982236324536
        0.213200720404175 -0.852802862171875 -0.426401436846280 -0.213200717039195
        -0.670820390731120 -0.670820393679672 -0.223606803382126 0.223606798385076
        0.684712473953190 -0.576599971377095 -0.432449984589644 -0.108112496267792
        0 0.725240663217528 0.644658376027243 0.241746889605056
        0.845154255153519 0.507092551317918 0.169030853378268 0
        -0.476731293667808 0.572077549477169 0.667423816642780 0
        0.966987556423525 -0.080582295755213 -0.241746891051114 0
        0.190692518611844 0.858116329163754 0.476731296400809 0
        -0.891882583595591 0.382235396183497 0.229341237473420 0.076447078279419
        0.452267019667424 0.753778359503679 0.452267016497876 0.150755674696506
        -0.725240663217528 -0.644658376027243 -0.241746889605056 0
        0.667423816642780 -0.572077549477169 -0.476731293667808 0
        0 -0.858116334300252 -0.476731287328518 -0.190692518178328
        0.169030857377403 -0.845154257382414 -0.507092546270047 0
        0.966987555641265 0.241746894314021 -0.080582295353615 0
        -0.753778361807378 -0.452267017803684 -0.452267014823439 -0.150755673792545
        0.632455530700779 -0.632455535202518 -0.316227759128976 -0.316227769232811
        -0.484164824942367 -0.789953151875243 -0.331270670560969 0.178376520195296
        0.801783728366898 -0.267261238306888 -0.534522481683180 -0.000000000141890
        0 -0.667423816642780 -0.572077549477169 -0.476731293667808
        0.241746889605056 -0.644658376027242 -0.725240663217528 0
        0.845154255153519 0.507092551317918 -0.169030853378267 0
        0 0.624695050559008 0.624695050645778 0.468521277537901
        0 -0.883452212909251 -0.331294568553883 -0.331294576400586
        0 0.780868811482643 0.624695045004908 0
        -0.468521279391265 0.624695044006864 0.624695055807899 0
        0 0.997054486958004 0.076696479951245 0
        0.970142498610814 0 0 -0.242535631174406
        0 0.948683302054730 0.316227754004188 0
        0 0.843661491133462 0.536875486848173 0
        0.997054484881728 0.076696506942838 0 0
        0.843661480036632 0.536875504286049 0 0
        0 -0.762492849732404 -0.457495709456026 -0.457495715757301
        0 -0.646996633565712 -0.539163877174602 -0.539163861645684
        0.762492854380984 -0.457495699618160 -0.457495717847533 0
        0 -0.833907847520107 -0.530668632008810 -0.151619605745038
        -0.557086012260807 -0.742781352512219 -0.371390680151605 0
        0.816496583028411 -0.408248288866001 -0.408248287860354 0
        0 0.742781353495836 0.557086015420803 0.371390673444375
        0.833907850067280 -0.151619607976834 -0.530668627368452 0
        0.408248288759059 0.816496580457949 0.408248293108221 0
        -0.964901283547497 0.214422499965823 0.107211252903483 0.107211248326215
        0.262612857660079 0.787838598999781 0.393919297747376 0.393919301100978
        0 0.974391195298725 0.224859508414302 0
        0.974391193141513 0.224859517762222 0 0
        -0.547722560173078 -0.730296742050804 -0.365148373038041 -0.182574180253128
        0.816496580386193 -0.408248293167179 -0.408248288843613 -0.000000002558343
        0 -0.745355997130286 -0.596284791467231 -0.298142390489507
        0.333333331564694 -0.666666672727882 -0.666666661489771 0
        0.894427193829204 0.447213589841382 0.000000000634355 0
        0.833333331383601 0.500000000937962 0.166666668616399 0.166666671651709
        -0.968962788295739 0.074535600304508 0.223606805515039 -0.074535600304508
        0.166666671651709 0.833333331383601 0.500000000937962 0.166666668616399
        -0.521749197090583 0.670820389263352 0.521749197090583 0.074535602360415
        -0.745355987202233 0.596284797395090 0.298142403453921 0
        0.666666672589671 0.666666662767202 0.333333329286253 0
        0 -0.963624112514936 -0.222374790221296 -0.148249864936997
        -0.832050295421840 -0.554700194599235 0.000000004677507 0
        0.534522481692390 -0.801783728187100 -0.267261238827861 0
        0.963624112514936 -0.222374790221296 -0.148249864936997 0
        0.267261238827861 0.801783728187100 0.534522481692390 0
        0.963624113112349 0.148249853042851 -0.222374795561934 0
        0.000000009286853 0.832050294181909 0.554700196459131 0
        0 -0.806559138596189 -0.513264900437760 -0.293294217340862
        -0.762000761773908 -0.635000633564652 -0.127000135545102 0
        0.577350275942079 -0.577350272669120 -0.577350258957679 0
        0 0.762000761773908 0.635000633564652 0.127000135545102
        0.806559138596189 0.513264900437760 0.293294217340862 0
        -0.577350258957679 0.577350272669120 0.577350275942079 0
        -0.933256525420056 0.207390343594001 0.207390333611325 0.207390338900900
        0.359210603630913 0.538815912686332 0.538815898655996 0.538815907180404
        0 -0.725476246947308 -0.652928629001294 -0.217642873868496
        0.229415733083201 -0.688247198540357 -0.688247204945467 0
        0.948683298508074 0.316227764644158 0.000000001851076 0
        0.974679434237846 -0.153896755372651 -0.153896751831812 -0.051298917485459
        -0.223606798809410 -0.670820389653621 -0.670820395626604 -0.223606800349493
        -0.725476252760983 0.652928621021009 0.217642878430436 0
        0.688247198817316 0.688247202472182 0.229415739672180 0
        0 0.933345605901459 0.358979080093030 0
        0.933345606668893 0.358979078097704 0 0
        0.534522479742128 -0.801783728459087 -0.267261241912424 -0.000000001360907
        0 -0.857142857142857 -0.428571428571429 -0.285714285714286
        0.303045756332547 -0.808122038767598 -0.505076271494648 0
        0.909137292441362 0.404061011569225 -0.101015258362508 0
        -0.707106786656749 -0.606091522042582 -0.303045761021291 -0.202030507347527
        -0.428571428571429 -0.857142857142857 -0.285714285714286 0.000000000000000
        0.857142857142857 -0.285714285714286 -0.428571428571429 0.000000000000000
        -0.801783725737273 -0.534522482463942 -0.267261244634238 0.000000004082720
        0 -0.808122038767598 -0.505076271494648 -0.303045756332547
        -0.285714285714286 -0.857142857142857 -0.428571428571428 0
        0.909137292441362 -0.101015258362508 -0.404061011569225 0
        0.944911183022734 -0.188982237004280 -0.188982236004947 -0.188982234006281
        -0.606091527512784 0.707106781186547 0.303045766491493 0.202030501877325
        0.785714285714286 0.500000000000000 0.357142857142857 0.071428571428572
        0.909137291659904 0.404061015476512 0.000000005470202 -0.101015249766477
        -0.357142857142857 0.785714285714286 0.500000000000000 -0.071428571428571
        0 0.909137287752617 0.404061024072544 0.101015250547934
        -0.857142857142857 0.428571428571429 0.285714285714286 0
        0.505076273057563 0.808122032515939 0.303045770398780 0
        -0.816496585467462 -0.408248286132412 -0.408248285715841 -0.000000001151214
        0.522232961316233 -0.696310626380971 -0.348155312513861 -0.348155316018862
        0 -0.923869770810735 -0.355334524810598 -0.142133817438871
        -0.301511345705917 -0.904534034783616 -0.301511340298633 0
        0.942809040236070 -0.235702260056406 -0.235702266118596 0
        -0.552770799765969 0.753778358660623 0.251259457755904 0.251259455203577
        0.833333332422321 0.500000004381119 0.166666664528566 0.166666660216471
        0.904534033174848 -0.402015126826162 -0.100503779984344 -0.100503785203577
        -0.426401433895857 -0.852802865517898 -0.213200718102617 -0.213200711857504
        -0.738548947783525 -0.615457451693498 -0.246182984646560 0.123091490174196
        0 0.773957301361971 0.633237787618913 0
        0.773957299669100 0.633237789687977 0 0
        0 0.700140048995795 0.700140035555833 0.140028005784712
        -0.408248289906120 -0.816496580354989 -0.408248292167081 -0.000000000295722
        0.857492929609518 -0.342997171707444 -0.171498584540379 -0.342997159421221
        0 -0.980196057998110 -0.140028020280892 -0.140028002275187
        0 0.693103277638181 0.693103281367391 0.198029512661537
        -0.408248292167081 -0.816496580354989 -0.408248289906120 -0.000000000295722
        0.848874685429528 -0.363803439081242 -0.121267811221093 -0.363803441587934
        0 -0.990147543346928 -0.099014749626712 -0.099014755266086
        0.921954446867191 -0.379628299078095 -0.054232613933750 -0.054232610363115
        -0.387298331911993 -0.903696116087519 -0.129099440491905 -0.129099443574320
        -0.708491907984497 -0.619930422274590 -0.309965206128506 0.132842232925558
        0.632455533790607 -0.632455532944699 -0.316227764678309 -0.316227762019457
        -0.387298333332206 -0.645497232456315 -0.645497217268248 -0.129099443795398
        0.921954446270582 -0.271163070429574 -0.271163071896235 -0.054232616356838
        0 0.920357988314051 0.383482490306811 0.076696499102633
        0.980196059194896 -0.140028006986959 -0.140028007191619 0
        0.182574183728686 0.912870931121890 0.365148367856762 0
        -0.532290646448842 0.720157934262850 0.344423360083532 0.281800934083103
        0.816496580579257 0.408248290518640 0.408248291106025 -0.000000003681816
        0.670820390988869 0.670820394026744 0.223606798612852 -0.223606801339887
        -0.685994341509776 0.514495753914047 0.514495755688018 0.000000003610203
        0 -0.690268486234828 -0.613571994821754 -0.383482495143563
        0.140028009084006 -0.700140037652055 -0.700140046239714 0
        0.912870928768376 0.365148373742829 -0.182574183724126 0
        -0.140028005784712 0.700140035555833 0.700140048995795 0
        0.980196056489980 -0.140028020554422 -0.140028012558565 0
        0.980196059522041 0.140028006831371 -0.140028005057196 0
        -0.076696497452289 0.920357984208100 0.383482500491163 0
        -0.939336434387938 0.313112151934966 0.125244859749666 0.062622427678073
        0.342997176418924 0.857492921742994 0.342997174128190 0.171498585036103
        -0.700140033125509 -0.700140049940932 -0.140028013210643 0
        0.690268498452352 -0.613571981081822 -0.383482495135911 0
        0 0.759072111297100 0.552052456556453 0.345032773893942
        0.845154256952159 -0.169030848149755 -0.507092550063023 0
        0.408248290743142 0.816496575365604 0.408248301308828 0
        -0.763762615533679 -0.545544725620123 -0.327326838496534 -0.109108937585709
        0.632455531364151 -0.632455531132883 -0.316227769781997 -0.316227765392314
        0 -0.897085228686456 -0.345032779691852 -0.276026218701562
        -0.507092551878789 -0.845154254729393 -0.169030853816286 0
        0.816496583225169 -0.408248290444576 -0.408248285888264 0
        -0.903696114731673 -0.387298333921916 0.129099437644571 0.129099449882804
        0.267261243641462 -0.801783723306901 -0.534522486605889 -0.000000002138902
        -0.878310066527975 0.390360030452349 0.195180013021381 0.195180009677830
        0.478091442127574 0.717137166078063 0.358568583212714 0.358568583574647
        -0.292770025240100 0.780720053832082 0.390360027913988 0.390360036981672
        0.956182886440095 0.239045723400798 0.119522859965539 0.119522867052862
        0 0.961523948380952 0.274721125306928 0
        0 0.880471106573440 0.474099810682669 0
        0.880471098763109 0.474099825187570 0 0
        0 0.738271659678350 0.671156056089019 0.067115608399302
        0.813733476752890 0.464990550375249 0.348742909423800 0
        -0.577350261038172 0.577350271877269 0.577350274653436 0
        0 -0.813733476752890 -0.464990550375249 -0.348742909423800
        -0.738271659678350 -0.671156056089019 -0.067115608399302 0
        0.577350274653436 -0.577350271877269 -0.577350261038172 0
        0.996615895645081 -0.047457899800570 -0.047457891406020 -0.047457905952176
        0.082199492380135 0.575396456035458 0.575396458037437 0.575396452815148];

    nA = [-0.948683300167071 -0.316227759667166 0.000000000000000
        -1.000000000000000 0 -0.000000000000000
        -0.577350267951537 -0.577350271803421 -0.577350267813920
        -0.816496581900497 -0.408248283238622 -0.408248295743562
        -0.707106779149207 -0.707106783223888 0.000000008792256
        -0.894427182880770 -0.447213611738249 0.000000000000000
        -0.980580674639535 -0.196116140395112 0
        -1.000000000000000 0 -0.000000000000000
        -0.912870925491948 -0.365148377716096 -0.182574192159729
        -0.912870925491948 -0.365148377716096 -0.182574192159729
        -0.912870925491948 -0.365148377716096 -0.182574192159729
        -0.894427191645603 -0.447213594208583 -0.000000003880661
        -0.894427191645603 -0.447213594208583 -0.000000003880661
        -0.816496584025179 -0.408248291221410 -0.408248283511409
        -0.816496584025179 -0.408248291221410 -0.408248283511409
        -0.857492918167444 -0.514495768002693 0.000000000000000
        -1.000000000000000 0 -0.000000000000000
        -0.666666675251644 -0.666666660164754 -0.333333329167204
        -0.666666675251644 -0.666666660164754 -0.333333329167204
        -0.942809043055013 -0.235702255875100 -0.235702259024135
        -0.942809043055013 -0.235702255875100 -0.235702259024135
        -0.666666675251644 -0.666666660164754 -0.333333329167204
        -0.942809043055013 -0.235702255875100 -0.235702259024135
        -0.707106786427146 -0.707106775945949 0.000000001667940
        -0.771516756494377 -0.617213391055954 -0.154303351712157
        -0.771516756494377 -0.617213391055954 -0.154303351712157
        -0.771516756494377 -0.617213391055954 -0.154303351712157
        -0.801783725639041 -0.534522487575102 -0.267261234706613
        -0.801783725639041 -0.534522487575102 -0.267261234706613
        -0.801783725639041 -0.534522487575102 -0.267261234706613
        -0.577350272548962 -0.577350264513355 -0.577350270506560
        -0.577350272548962 -0.577350264513355 -0.577350270506560
        -0.904534031595092 -0.301511355883906 -0.301511339686217
        -0.639602149757223 -0.639602145352570 -0.426401437247025
        -0.707106779725268 -0.707106782647827 -0.000000011688130
        -0.989949493510053 -0.141421357295107 0
        -0.989949493510053 -0.141421357295107 0
        -0.989949493510053 -0.141421357295107 0
        -1.000000000000000 0 0.000000000000000
        -0.832050294073944 -0.554700196621079 -0.000000000000000
        -0.919145032908829 -0.393919291834036 0
        -1.000000000000000 0 -0.000000000000000
        -0.861640435364518 -0.492365966086633 -0.123091492737805
        -0.861640435364518 -0.492365966086633 -0.123091492737805
        -0.861640435364518 -0.492365966086633 -0.123091492737805
        -0.904534028759875 -0.301511352534824 -0.301511351540951
        -0.904534028759875 -0.301511352534824 -0.301511351540951
        -0.816496583953726 -0.408248285838872 -0.408248289036853
        -0.816496583953726 -0.408248285838872 -0.408248289036853
        -0.727606878247232 -0.485071240616695 -0.485071254821287
        -0.685994337241436 -0.514495755088474 -0.514495760204711
        -0.970142502121177 -0.242535617132954 -0.000000000000000
        -0.707106784197378 -0.707106778175717 0.000000009789942
        -0.717137165764438 -0.597614307820115 -0.358568577217851
        -0.717137165764438 -0.597614307820115 -0.358568577217851
        -0.717137165764438 -0.597614307820115 -0.358568577217851
        -0.894427194111714 -0.447213589276361 -0.000000006386481
        -0.894427194111714 -0.447213589276361 -0.000000006386481
        -0.801783731495532 -0.534522474911388 -0.267261242464568
        -0.801783731495532 -0.534522474911388 -0.267261242464568
        -0.813733472837210 -0.581238191436432 0
        -1.000000000000000 0 -0.000000000000000
        -0.688247204012338 -0.688247198794377 -0.229415735120529
        -0.973328526489955 -0.162221427587113 -0.162221416442133
        -0.707106788035655 -0.707106774337441 0.000000005970075
        -0.792593921867616 -0.566138517742786 -0.226455412270304
        -0.792593921867616 -0.566138517742786 -0.226455412270304
        -0.792593921867616 -0.566138517742786 -0.226455412270304
        -0.784464543605915 -0.588348403541864 -0.196116128543533
        -0.784464543605915 -0.588348403541864 -0.196116128543533
        -0.784464543605915 -0.588348403541864 -0.196116128543533
        -0.577350269615320 -0.577350264539909 -0.577350273413649
        -0.577350269615320 -0.577350264539909 -0.577350273413649
        -0.993883735904950 -0.110431514992866 0.000000000000000
        -1.000000000000000 0 -0.000000000000000
        -0.872871563406380 -0.436435781562180 -0.218217878205962
        -0.872871563406380 -0.436435781562180 -0.218217878205962
        -0.872871563406380 -0.436435781562180 -0.218217878205962
        -0.801783728406615 -0.534522485851651 -0.267261229850796
        -0.801783728406615 -0.534522485851651 -0.267261229850796
        -0.816496582303071 -0.408248294066437 -0.408248284110599
        -0.816496582303071 -0.408248294066437 -0.408248284110599
        -0.737864781795359 -0.527046281046556 -0.421637025676211
        -0.737864781795359 -0.527046281046556 -0.421637025676211
        -0.737864781795359 -0.527046281046556 -0.421637025676211
        -0.737864781795359 -0.527046281046556 -0.421637025676211
        -0.737864781795359 -0.527046281046556 -0.421637025676211
        -0.843274042330571 -0.527046276844005 -0.105409257640543
        -0.843274042330571 -0.527046276844005 -0.105409257640543
        -0.843274042330571 -0.527046276844005 -0.105409257640543
        -0.843274042330571 -0.527046276844005 -0.105409257640543
        -0.843274042330571 -0.527046276844005 -0.105409257640543
        -0.843274042330571 -0.527046276844005 -0.105409257640543
        -0.843274042330571 -0.527046276844005 -0.105409257640543
        -0.843274042330571 -0.527046276844005 -0.105409257640543
        -0.843274042330571 -0.527046276844005 -0.105409257640543
        -0.666666666666667 -0.666666666666667 -0.333333333333333
        -0.666666666666667 -0.666666666666667 -0.333333333333333
        -0.666666666666667 -0.666666666666667 -0.333333333333333
        -0.666666666666667 -0.666666666666667 -0.333333333333333
        -0.666666666666667 -0.666666666666667 -0.333333333333333
        -0.800000007589466 -0.599999989880711 0.000000000000000
        -0.800000007589466 -0.599999989880711 0.000000000000000
        -0.800000002904363 -0.599999996127516 -0.000000011686254
        -0.800000002904363 -0.599999996127516 -0.000000011686254
        -0.989949492443151 -0.141421364763421 0
        -0.989949492443151 -0.141421364763421 0
        -0.707106778519293 -0.565685424178679 -0.424264074184765
        -0.707106778519293 -0.565685424178679 -0.424264074184765
        -0.707106778519293 -0.565685424178679 -0.424264074184765
        -0.707106778519293 -0.565685424178679 -0.424264074184765
        -0.874157276519123 -0.485642930462977 0
        -1.000000000000000 0 0.000000000000000
        -0.962250448705994 -0.192450081347409 -0.192450097829253
        -0.962250448705994 -0.192450081347409 -0.192450097829253
        -0.962250448705994 -0.192450081347409 -0.192450097829253
        -0.680413820500869 -0.680413814458744 -0.272165526775735
        -0.680413820500869 -0.680413814458744 -0.272165526775735
        -0.680413820500869 -0.680413814458744 -0.272165526775735
        -0.962250450855370 -0.192450082847319 -0.192450085582461
        -0.962250450855370 -0.192450082847319 -0.192450085582461
        -0.952579344260516 -0.272165528969358 -0.136082760587207
        -0.952579344260516 -0.272165528969358 -0.136082760587207
        -0.952579344260516 -0.272165528969358 -0.136082760587207
        -0.952579344260516 -0.272165528969358 -0.136082760587207
        -0.952579344260516 -0.272165528969358 -0.136082760587207
        -0.952579344260516 -0.272165528969358 -0.136082760587207
        -0.952579344260516 -0.272165528969358 -0.136082760587207
        -0.942809045169483 -0.235702252396786 -0.235702254044567
        -0.942809045169483 -0.235702252396786 -0.235702254044567
        -0.942809045169483 -0.235702252396786 -0.235702254044567
        -0.707106781276638 -0.707106781096458 -0.000000002146545
        -0.953462587590819 -0.286038780861493 -0.095346263208841
        -0.953462587590819 -0.286038780861493 -0.095346263208841
        -0.953462587590819 -0.286038780861493 -0.095346263208841
        -0.904534031840555 -0.301511351744986 -0.301511343088750
        -0.904534031840555 -0.301511351744986 -0.301511343088750
        -0.948683299991224 -0.316227760194705 -0.000000008875729
        -0.948683299991224 -0.316227760194705 -0.000000008875729
        -0.749268645059556 -0.655610072251617 -0.093658585793540
        -0.749268645059556 -0.655610072251617 -0.093658585793540
        -0.749268645059556 -0.655610072251617 -0.093658585793540
        -0.811107102130947 -0.486664267127736 -0.324442845465517
        -0.811107102130947 -0.486664267127736 -0.324442845465517
        -0.811107102130947 -0.486664267127736 -0.324442845465517
        -0.577350270717785 -0.577350279502621 -0.577350257348472
        -0.577350270717785 -0.577350279502621 -0.577350257348472
        -0.928476690782124 -0.371390676611943 -0.000000000000000
        -0.995893206728755 -0.090535743170954 -0.000000000000000
        -1.000000000000000 0 -0.000000000000000
        -0.964763821767187 -0.263117403850776 -0.000000000000000
        -0.964763821767187 -0.263117403850776 -0.000000000000000
        -0.964763821767187 -0.263117403850776 -0.000000000000000
        -0.789352223846899 -0.613940605195613 -0.000000000000000
        -1.000000000000000 0 -0.000000000000000
        -1.000000000000000 0 0.000000000000000
        -0.696310618597245 -0.696310628500968 -0.174077658145175
        -0.870388277890326 -0.348155319575241 -0.348155308967958
        -0.870388277890326 -0.348155319575241 -0.348155308967958
        -0.870388277890326 -0.348155319575241 -0.348155308967958
        -0.984731927447613 -0.123091491129882 -0.123091493925171
        -0.984731927447613 -0.123091491129882 -0.123091493925171
        -0.984731927447613 -0.123091491129882 -0.123091493925171
        -0.615457452585078 -0.615457454537983 -0.492365967255163
        -0.707106785819505 -0.707106776553590 0.000000003165454
        -0.707106783777999 -0.707106778595097 0.000000005409034
        -0.845154252121480 -0.507092559738786 -0.169030843275855
        -0.845154252121480 -0.507092559738786 -0.169030843275855
        -0.845154252121480 -0.507092559738786 -0.169030843275855
        -0.801783727099916 -0.534522478536335 -0.267261248401523
        -0.801783727099916 -0.534522478536335 -0.267261248401523
        -0.948683297520349 -0.316227767607332 -0.000000007829098
        -0.948683297520349 -0.316227767607332 -0.000000007829098
        -0.845154252121480 -0.507092559738786 -0.169030843275855
        -0.845154252121480 -0.507092559738786 -0.169030843275855
        -0.845154252121480 -0.507092559738786 -0.169030843275855
        -0.845154252121480 -0.507092559738786 -0.169030843275855
        -0.845154252121480 -0.507092559738786 -0.169030843275855
        -0.845154252121480 -0.507092559738786 -0.169030843275855
        -0.845154259642273 -0.507092540914963 -0.169030862143363
        -0.845154259642273 -0.507092540914963 -0.169030862143363
        -0.845154259642273 -0.507092540914963 -0.169030862143363
        -0.845154259642273 -0.507092540914963 -0.169030862143363
        -0.845154255307662 -0.507092553773183 -0.169030845241758
        -0.845154255307662 -0.507092553773183 -0.169030845241758
        -0.979957885397199 -0.178174163495740 -0.089087093961156
        -0.979957885397199 -0.178174163495740 -0.089087093961156
        -0.979957885397199 -0.178174163495740 -0.089087093961156
        -0.979957885397199 -0.178174163495740 -0.089087093961156
        -0.979957885397199 -0.178174163495740 -0.089087093961156
        -0.979957885397199 -0.178174163495740 -0.089087093961156
        -0.979957885397199 -0.178174163495740 -0.089087093961156
        -0.979957885397199 -0.178174163495740 -0.089087093961156
        -0.979957885397199 -0.178174163495740 -0.089087093961156
        -0.979957885397199 -0.178174163495740 -0.089087093961156
        -0.979957885397199 -0.178174163495740 -0.089087093961156
        -0.890870803965847 -0.445435406252094 -0.089087089402880
        -0.890870803965847 -0.445435406252094 -0.089087089402880
        -0.890870803965847 -0.445435406252094 -0.089087089402880
        -0.890870803965847 -0.445435406252094 -0.089087089402880
        -0.890870803965847 -0.445435406252094 -0.089087089402880
        -0.890870803965847 -0.445435406252094 -0.089087089402880
        -0.890870803965847 -0.445435406252094 -0.089087089402880
        -0.843274043185669 -0.527046275291862 -0.105409258560479
        -0.843274043185669 -0.527046275291862 -0.105409258560479
        -0.843274043185669 -0.527046275291862 -0.105409258560479
        -0.843274043185669 -0.527046275291862 -0.105409258560479
        -0.843274043185669 -0.527046275291862 -0.105409258560479
        -0.843274043185669 -0.527046275291862 -0.105409258560479
        -0.737864790476569 -0.527046274415987 -0.421637018772305
        -0.737864790476569 -0.527046274415987 -0.421637018772305
        -0.737864790476569 -0.527046274415987 -0.421637018772305
        -0.737864790476569 -0.527046274415987 -0.421637018772305
        -0.910366473688846 -0.413802952603453 0
        -1.000000000000000 0 -0.000000000000000
        -0.986393924174848 -0.164398985249130 0.000000000000000
        -0.966987555641265 -0.241746894314021 -0.080582295353615
        -0.966987555641265 -0.241746894314021 -0.080582295353615
        -0.966987555641265 -0.241746894314021 -0.080582295353615
        -0.966987555641265 -0.241746894314021 -0.080582295353615
        -0.966987555641265 -0.241746894314021 -0.080582295353615
        -0.966987555641265 -0.241746894314021 -0.080582295353615
        -0.966987555641265 -0.241746894314021 -0.080582295353615
        -0.725240663217528 -0.644658376027243 -0.241746889605056
        -0.725240663217528 -0.644658376027243 -0.241746889605056
        -0.725240663217528 -0.644658376027243 -0.241746889605056
        -0.845154253198270 -0.507092554898358 -0.169030852413195
        -0.845154253198270 -0.507092554898358 -0.169030852413195
        -0.845154253198270 -0.507092554898358 -0.169030852413195
        -0.845154253198270 -0.507092554898358 -0.169030852413195
        -0.845154255153519 -0.507092551317918 -0.169030853378268
        -0.845154255153519 -0.507092551317918 -0.169030853378268
        -0.858116334300252 -0.476731287328518 -0.190692518178328
        -0.858116334300252 -0.476731287328518 -0.190692518178328
        -0.858116334300252 -0.476731287328518 -0.190692518178328
        -0.858116334300252 -0.476731287328518 -0.190692518178328
        -0.858116334300252 -0.476731287328518 -0.190692518178328
        -0.858116334300252 -0.476731287328518 -0.190692518178328
        -0.858116334300252 -0.476731287328518 -0.190692518178328
        -0.667423816642780 -0.572077549477169 -0.476731293667808
        -0.667423816642780 -0.572077549477169 -0.476731293667808
        -0.667423816642780 -0.572077549477169 -0.476731293667808
        -0.624695050559008 -0.624695050645778 -0.468521277537901
        -0.883452212909251 -0.331294568553883 -0.331294576400586
        -0.780868811482643 -0.624695045004908 -0.000000000000000
        -0.707106777903288 -0.707106784469807 -0.000000009055127
        -0.997054486958004 -0.076696479951245 0.000000000000000
        -0.997054486958004 -0.076696479951245 0.000000000000000
        -0.997054486958004 -0.076696479951245 0.000000000000000
        -0.843661491133462 -0.536875486848173 0.000000000000000
        -1.000000000000000 0 0.000000000000000
        -1.000000000000000 0 0.000000000000000
        -0.762492849732404 -0.457495709456026 -0.457495715757302
        -0.646996633565712 -0.539163877174602 -0.539163861645684
        -0.707106777024286 -0.707106785348809 -0.000000011910511
        -0.833907847520107 -0.530668632008810 -0.151619605745039
        -0.833907847520107 -0.530668632008810 -0.151619605745039
        -0.833907847520107 -0.530668632008810 -0.151619605745039
        -0.742781353495836 -0.557086015420803 -0.371390673444375
        -0.742781353495836 -0.557086015420803 -0.371390673444375
        -0.742781353495836 -0.557086015420803 -0.371390673444375
        -0.816496578703176 -0.408248291157391 -0.408248294219436
        -0.816496578703176 -0.408248291157391 -0.408248294219436
        -0.974391195298725 -0.224859508414302 -0.000000000000000
        -1.000000000000000 0 0.000000000000000
        -0.745355997130286 -0.596284791467231 -0.298142390489507
        -0.745355997130286 -0.596284791467231 -0.298142390489507
        -0.745355997130286 -0.596284791467231 -0.298142390489507
        -0.745355997130286 -0.596284791467231 -0.298142390489507
        -0.745355997130286 -0.596284791467231 -0.298142390489507
        -0.745355987202233 -0.596284797395090 -0.298142403453921
        -0.745355987202233 -0.596284797395090 -0.298142403453921
        -0.745355987202233 -0.596284797395090 -0.298142403453921
        -0.745355987202233 -0.596284797395090 -0.298142403453921
        -0.894427190356857 -0.447213596786075 -0.000000004421533
        -0.894427190356857 -0.447213596786075 -0.000000004421533
        -0.963624112514936 -0.222374790221296 -0.148249864936997
        -0.963624112514936 -0.222374790221296 -0.148249864936997
        -0.963624112514936 -0.222374790221296 -0.148249864936997
        -0.832050295421840 -0.554700194599235 -0.000000004677507
        -0.832050295421840 -0.554700194599235 -0.000000004677507
        -0.801783725419136 -0.534522486921338 -0.267261236673856
        -0.801783725419136 -0.534522486921338 -0.267261236673856
        -0.806559138596189 -0.513264900437760 -0.293294217340862
        -0.806559138596189 -0.513264900437760 -0.293294217340862
        -0.806559138596189 -0.513264900437760 -0.293294217340862
        -0.762000761773908 -0.635000633564652 -0.127000135545102
        -0.762000761773908 -0.635000633564652 -0.127000135545102
        -0.762000761773908 -0.635000633564652 -0.127000135545102
        -0.577350278106534 -0.577350261731361 -0.577350267730983
        -0.577350278106534 -0.577350261731361 -0.577350267730983
        -0.725476246947308 -0.652928629001294 -0.217642873868497
        -0.725476246947308 -0.652928629001294 -0.217642873868497
        -0.725476246947308 -0.652928629001294 -0.217642873868497
        -0.688247198540357 -0.688247204945467 -0.229415733083201
        -0.688247198540357 -0.688247204945467 -0.229415733083201
        -0.948683295868110 -0.316227772564048 0.000000000326930
        -0.948683295868110 -0.316227772564048 0.000000000326930
        -0.933345605901459 -0.358979080093030 0.000000000000000
        -1.000000000000000 0 0.000000000000000
        -0.857142857142857 -0.428571428571429 -0.285714285714286
        -0.857142857142857 -0.428571428571429 -0.285714285714286
        -0.857142857142857 -0.428571428571429 -0.285714285714286
        -0.857142857142857 -0.428571428571429 -0.285714285714286
        -0.857142857142857 -0.428571428571429 -0.285714285714286
        -0.857142857142857 -0.428571428571429 -0.285714285714286
        -0.857142857142857 -0.428571428571429 -0.285714285714286
        -0.808122038767598 -0.505076271494648 -0.303045756332547
        -0.808122038767598 -0.505076271494648 -0.303045756332547
        -0.808122038767598 -0.505076271494648 -0.303045756332547
        -0.808122038767598 -0.505076271494648 -0.303045756332547
        -0.909137287752617 -0.404061024072544 -0.101015250547934
        -0.909137287752617 -0.404061024072544 -0.101015250547934
        -0.909137287752617 -0.404061024072544 -0.101015250547934
        -0.909137287752617 -0.404061024072544 -0.101015250547934
        -0.909137287752617 -0.404061024072544 -0.101015250547934
        -0.909137287752617 -0.404061024072544 -0.101015250547934
        -0.909137287752617 -0.404061024072544 -0.101015250547934
        -0.909137287752617 -0.404061024072544 -0.101015250547934
        -0.923869770810735 -0.355334524810598 -0.142133817438871
        -0.923869770810735 -0.355334524810598 -0.142133817438871
        -0.923869770810735 -0.355334524810598 -0.142133817438871
        -0.923869770810735 -0.355334524810598 -0.142133817438871
        -0.923869770810735 -0.355334524810598 -0.142133817438871
        -0.904534034493290 -0.301511339558968 -0.301511347316561
        -0.904534034493290 -0.301511339558968 -0.301511347316561
        -0.942809043716874 -0.235702256518229 -0.235702255733560
        -0.942809043716874 -0.235702256518229 -0.235702255733560
        -0.942809043716874 -0.235702256518229 -0.235702255733560
        -0.773957301361971 -0.633237787618913 -0.000000000000000
        -1.000000000000000 0 0.000000000000000
        -0.700140048995795 -0.700140035555833 -0.140028005784712
        -0.980196057998110 -0.140028020280892 -0.140028002275187
        -0.980196057998110 -0.140028020280892 -0.140028002275187
        -0.980196057998110 -0.140028020280892 -0.140028002275187
        -0.693103277638181 -0.693103281367391 -0.198029512661537
        -0.990147543346928 -0.099014749626712 -0.099014755266086
        -0.990147543346928 -0.099014749626712 -0.099014755266086
        -0.990147543346928 -0.099014749626712 -0.099014755266086
        -0.980196059194896 -0.140028006986959 -0.140028007191619
        -0.980196059194896 -0.140028006986959 -0.140028007191619
        -0.980196059194896 -0.140028006986959 -0.140028007191619
        -0.980196059194896 -0.140028006986959 -0.140028007191619
        -0.700140049940932 -0.700140033125509 -0.140028013210643
        -0.700140049940932 -0.700140033125509 -0.140028013210643
        -0.920357988314051 -0.383482490306811 -0.076696499102633
        -0.920357988314051 -0.383482490306811 -0.076696499102633
        -0.920357988314051 -0.383482490306811 -0.076696499102633
        -0.920357988314051 -0.383482490306811 -0.076696499102633
        -0.920357988314051 -0.383482490306811 -0.076696499102633
        -0.920357988314051 -0.383482490306811 -0.076696499102633
        -0.920357988314051 -0.383482490306811 -0.076696499102633
        -0.690268486234828 -0.613571994821754 -0.383482495143563
        -0.690268486234828 -0.613571994821754 -0.383482495143563
        -0.690268486234828 -0.613571994821754 -0.383482495143563
        -0.707106774137940 -0.707106788235155 0.000000002617597
        -0.707106778876181 -0.707106783496913 0.000000006428261
        -0.912870927139619 -0.365148378343979 -0.182574182665610
        -0.912870927139619 -0.365148378343979 -0.182574182665610
        -0.912870927139619 -0.365148378343979 -0.182574182665610
        -0.912870927139619 -0.365148378343979 -0.182574182665610
        -0.912870928138593 -0.365148373267686 -0.182574187823324
        -0.912870928138593 -0.365148373267686 -0.182574187823324
        -0.759072111297100 -0.552052456556453 -0.345032773893942
        -0.759072111297100 -0.552052456556453 -0.345032773893942
        -0.759072111297100 -0.552052456556453 -0.345032773893942
        -0.897085228686456 -0.345032779691852 -0.276026218701562
        -0.897085228686456 -0.345032779691852 -0.276026218701562
        -0.897085228686456 -0.345032779691852 -0.276026218701562
        -0.897085228686456 -0.345032779691852 -0.276026218701562
        -0.897085228686456 -0.345032779691852 -0.276026218701562
        -0.897085228686456 -0.345032779691852 -0.276026218701562
        -0.897085228686456 -0.345032779691852 -0.276026218701562
        -0.816496580686100 -0.408248294264277 -0.408248287146700
        -0.816496580686100 -0.408248294264277 -0.408248287146700
        -0.816496575365604 -0.408248290743142 -0.408248301308827
        -0.816496575365604 -0.408248290743142 -0.408248301308827
        -0.961523948380952 -0.274721125306928 0
        -0.880471106573440 -0.474099810682669 -0.000000000000000
        -1.000000000000000 0 -0.000000000000000
        -0.738271659678350 -0.671156056089019 -0.067115608399302
        -0.738271659678350 -0.671156056089019 -0.067115608399302
        -0.738271659678350 -0.671156056089019 -0.067115608399302
        -0.813733476752890 -0.464990550375249 -0.348742909423800
        -0.813733476752890 -0.464990550375249 -0.348742909423800
        -0.813733476752890 -0.464990550375249 -0.348742909423800
        -0.577350278273464 -0.577350267268553 -0.577350262026861
        -0.577350278273464 -0.577350267268553 -0.577350262026861];

    y = [1.2180
        0.9580
        0.0640
        0.8790
        1.0050
        1.2850
        1.2010
        0.8840
        1.1440
        0.9780
        0.8090
        1.2980
        0.6920
        1.2650
        1.2830
        1.1800
        0.9380
        0.8100
        1.1600
        0.6150
        0.8950
        1.1540
        1.2340
        1.3280
        1.0750
        1.2050
        1.0000
        1.1870
        1.1600
        0.9780
        0.4890
        0.4890
        0.4240
        0.9800
        1.2420
        1.0990
        1.3060
        1.3650
        0.8020
        1.0910
        1.3450
        0.9820
        1.1790
        1.2520
        0.6030
        1.1330
        0.8740
        1.3750
        1.1800
        0.9190
        0.8490
        1.2410
        0.9100
        1.2630
        1.1790
        0.6820
        0.9550
        1.1970
        1.2540
        1.0650
        1.0190
        0.8440
        1.1390
        0.9940
        1.3000
        1.2170
        1.2290
        0.9880
        1.1570
        1.2330
        0.9940
        0.4980
        0.4980
        1.0010
        0.7310
        1.0850
        0.8240
        0.9350
        1.2390
        0.4730
        1.1020
        1.0970
        1.2800
        0.9410
        1.1260
        1.1090
        1.1560
        1.2170
        0.7350
        1.2450
        1.2740
        1.1530
        1.1560
        1.1330
        1.0330
        1.1510
        1.2340
        1.2510
        1.1060
        1.1420
        1.1620
        1.2770
        0.9390
        1.2420
        0.9980
        1.2380
        1.0350
        1.0440
        1.0960
        0.9160
        1.0980
        1.2550
        0.9620
        0.9410
        0.9380
        0.9030
        1.1910
        1.0820
        0.4100
        1.2800
        1.3430
        1.3100
        1.1120
        1.2070
        1.3090
        1.3320
        0.9000
        1.2820
        1.3550
        1.3480
        1.3250
        1.2960
        1.2560
        1.2570
        0.9670
        1.2920
        1.2110
        1.2730
        1.2630
        0.8580
        1.1830
        1.0040
        1.1810
        0.9670
        0.9500
        0.4480
        0.4480
        1.3360
        0.9230
        0.6720
        1.2750
        1.2050
        1.3220
        0.8830
        0.9550
        0.7660
        1.0880
        0.9450
        0.4920
        0.7790
        0.9800
        0.8540
        0.8760
        0.7310
        1.1990
        1.1880
        1.0290
        0.9440
        0.9080
        1.1670
        1.1300
        1.3720
        1.0070
        1.1080
        1.2250
        0.9670
        0.9600
        0.5670
        1.2350
        1.1790
        1.3250
        1.1000
        1.1830
        1.2400
        1.1870
        1.3320
        1.3240
        1.2830
        0.8240
        1.1710
        1.3520
        1.2170
        1.3140
        1.0200
        1.0900
        1.2650
        1.2040
        1.0670
        1.1510
        0.9560
        1.2830
        1.3750
        0.9590
        1.2860
        1.2390
        1.1430
        1.1980
        1.2430
        1.3340
        1.1560
        1.2910
        1.3180
        1.0430
        1.3530
        0.9820
        1.1520
        1.2460
        1.2580
        1.3850
        1.0450
        1.1760
        1.2670
        1.0250
        1.2610
        1.3160
        1.2310
        1.2790
        1.2680
        1.3150
        1.0750
        1.2600
        0.9400
        1.1260
        1.3730
        1.2380
        1.0730
        1.1770
        1.2030
        1.1570
        0.9220
        1.1070
        1.0550
        0.8380
        0.7090
        0.8240
        1.2180
        0.8480
        1.2920
        1.2890
        1.1520
        0.6220
        0.9140
        0.9430
        0.7090
        0.9970
        1.0690
        0.9290
        0.9670
        1.2430
        1.2740
        0.3510
        1.1120
        1.2360
        1.2500
        0.9220
        1.2130
        0.9660
        1.3140
        1.1680
        1.0880
        1.2420
        1.1940
        1.1970
        1.1730
        1.2030
        1.1100
        1.1810
        0.9070
        1.0880
        1.2720
        0.8890
        1.2940
        1.1440
        1.2270
        1.1070
        0.9670
        0.9900
        1.2410
        1.0030
        0.4790
        0.4780
        1.2180
        1.2680
        1.0360
        1.2130
        1.0620
        1.1820
        1.1620
        1.3310
        0.9830
        1.0880
        1.0530
        1.0670
        1.2210
        1.3480
        1.0970
        1.2810
        1.2170
        1.2170
        1.0840
        1.3410
        0.9100
        1.3590
        1.2280
        1.2890
        1.1070
        1.2910
        1.3900
        0.9190
        0.9500
        1.3430
        1.2440
        1.2470
        1.3370
        1.3170
        1.3060
        1.3310
        1.2660
        1.3070
        0.7820
        0.7010
        0.9960
        0.7040
        0.8420
        0.9820
        1.1540
        0.7400
        0.9400
        0.9070
        1.2920
        1.3720
        1.3010
        1.2930
        1.2390
        1.3030
        1.3270
        1.2270
        1.3250
        1.3600
        0.9600
        1.2120
        0.9450
        1.1850
        1.0380
        1.1490
        1.1180
        1.2690
        1.2820
        1.4140
        1.3070
        1.1120
        1.0340
        0.9720
        1.2890
        1.2720
        0.3840
        1.2190
        1.2480
        0.8360
        1.1240
        0.8950
        1.2480
        1.0570
        1.1520
        1.3570
        1.2320
        1.1250
        1.2740
        1.2790
        0.9680
        0.7220
        1.2020
        1.0050
        1.1500
        0.8150
        0.9330
        0.4040
        0.4040
        ];
end

function o = om2oct(omA, omB, mA, epsijk)

    arguments
        omA(3, 3) double {mustBeFinite, mustBeReal}
        omB(3, 3) double {mustBeFinite, mustBeReal}
        mA(:, 3) double {mustBeFinite, mustBeReal} = [0 0 1]
        epsijk(1, 1) double = 1
    end

    % OM2OCT Convert orientation matrices (omA,omB) to octonions (o).
    % Orientation matrices of grain A and grain B and sample frame boundary
    % plane normal pointing outward from grain A towards grain B to octonion
    % with BP normal = [0 0 1]; epsijk == 1 (active rotation matrix), epsijk ==
    % -1 (passive rotation matrix)
    %--------------------------------------------------------------------------
    % Author(s): Sterling Baird
    % Date: 2020-12-05
    %
    % Inputs:
    %  omA, omB - Orientation matrices of grains A and B in sample reference
    %  frame, resp.
    %
    % Outputs:
    %   o - octonion, with BP normal = [0 0 1]
    %
    % Usage:
    %  o = om2oct(omA,omB)
    %
    % Dependencies:
    %  om2qu.m
    %  qmA2oct.m
    %
    % References:
    %  (1) Francis, T.; Chesser, I.; Singh, S.; Holm, E. A.; De Graef, M. A
    %  Geodesic Octonion Metric for Grain Boundaries. Acta Materialia 2019,
    %  166, 135–147. https://doi.org/10.1016/j.actamat.2018.12.034.

    %--------------------------------------------------------------------------
    % % https://www.mathworks.com/matlabcentral/answers/361219-convert-a-3d-array-to-a-2d-cell-array
    % omA = num2cell(omA, 3);
    % omB = num2cell(omB, 3);

    pA = om2qu(omA, epsijk);
    pB = om2qu(omB, epsijk);

    o = qmA2oct(pA, pB, mA, epsijk);

end

function [qm, nA] = oct2five(o, epsijk)

    arguments
        o(:, 8) double
        epsijk(1, 1) double = 1
    end

    % OCT2FIVE Convert octonions (o) to 5DOF coordinates (qm,nA).
    % octonion with BP normal = [0 0 1]. Misorientation quaternion and grain A
    % crystal frame boundary plane normal pointing outward from grain A towards grain B.
    %--------------------------------------------------------------------------
    % Author(s): Sterling Baird
    % Date: 2020-01-27
    %
    % Inputs:
    %   o - octonion, with BP normal = [0 0 1]
    %
    % Outputs:
    %  qm - Misorientation quaternion
    %  nA - boundary plane normal (grain A crystal frame) pointing from grain
    %  A to grain B
    %
    % Usage:
    %  [qm,nA] = oct2five(o)
    %
    % Dependencies:
    %  qinv.m
    %  qmA2five.m
    %
    % References:
    %  (1) Francis, T.; Chesser, I.; Singh, S.; Holm, E. A.; De Graef, M. A
    %  Geodesic Octonion Metric for Grain Boundaries. Acta Materialia 2019,
    %  166, 135–147. https://doi.org/10.1016/j.actamat.2018.12.034.
    %--------------------------------------------------------------------------
    npts = size(o, 1);
    qA = normr(o(:, 1:4));
    qB = normr(o(:, 5:8));
    mA = repmat([0 0 1], npts, 1);
    [qm, nA] = qmA2five(qA, qB, mA, epsijk);
end

function [qm, nA] = qmA2five(qA, qB, mA, epsijk)

    arguments
        qA(:, 4) double {mustBeFinite, mustBeReal}
        qB(:, 4) double {mustBeFinite, mustBeReal}
        mA(:, 3) double {mustBeFinite, mustBeReal}
        epsijk(1, 1) double = 1
    end

    % EUMA2FIVE Convert sample frame euler angles of grain A and grain B and
    % sample frame boundary plane normal pointing outward from grain A towards
    % grain B to misorientation quaternion and boundary plane normal (crystal
    % reference frame of grain A).
    %--------------------------------------------------------------------------
    % Author(s): Sterling Baird
    % Date: 2020-08-22
    %
    % Inputs:
    %  eA, eB - Euler angles of grains A and B in sample reference frame,
    %    resp.
    %  mA - boundary plane normal (sample reference frame) pointing from grain
    %  A to grain B
    %
    % Outputs:
    %  five -struct containing qm and nA
    %   --qm - misorientation quaternion (qB qA*)
    %
    %   --nA - boundary plane normal in grain A crystal frame vec[qA mA qA*]
    %
    % Usage:
    %  five = eumA2qnA(euA,euB,mA)
    %
    % Dependencies:
    %  qmult.m
    %  qinv_francis.m
    %  eu2qu.m
    %
    % References:
    %  [1] supplementary material of DOI: 10.1016/j.actamat.2018.12.034
    %
    %  [2]
    %  Seita, M.; Volpi, M.; Patala, S.; McCue, I.; Schuh, C. A.; Diamanti, M.
    %  V.; Erlebacher, J.; Demkowicz, M. J. A High-Throughput Technique for
    %  Determining Grain Boundary Character Non-Destructively in
    %  Microstructures with through-Thickness Grains. Npj Computational
    %  Materials 2016, 2, 16016.

    %--------------------------------------------------------------------------
    npts = size(qA, 1);
    %misorientation quaternion
    qm = qlab2qm(qA, qB, epsijk);

    %boundary plane normal
    [nA, nB] = deal(zeros(npts, 3));

    for i = 1:npts
        mAtmp = mA(i, :);
        qAtmp = qA(i, :);
        nA(i, :) = Lpr(qinv(qAtmp), mAtmp, epsijk); %sample frame to crystal frame?
    end

end

function qm = qlab2qm(qAlab, qBlab, epsijk)

    arguments
        qAlab(:, 4) double
        qBlab(:, 4) double
        epsijk(1, 1) double = 1
    end

    % QLAB2QM  Convert lab/sample frame quaternions of grain A and grain B to misorientation quaternion
    % active (epsijk==1) or passive (epsijk==-1) convention
    %--------------------------------------------------------------------------
    % Author(s): Sterling Baird
    % Date: 2020-08-22
    %
    % Inputs:
    %  qA, qB - Orientations of grains A and B in sample reference frame, resp.
    %  convention - francis or johnson convention
    %
    % Outputs:
    %  qm - misorientation quaternion (francis: qB qA*, johnson: qA* qB)
    %
    % Usage:
    %  qm = qlab2qm(qAlab,qBlab,convention)
    %
    % Dependencies:
    %  qmult.m
    %  qinv.m
    %
    % References:
    %  [1] supplementary material of DOI: 10.1016/j.actamat.2018.12.034
    %
    %  [2] Rowenhorst, D.; Rollett, A. D.; Rohrer, G. S.; Groeber, M.; Jackson,
    %  M.; Konijnenberg, P. J.; De Graef, M. Consistent Representations of and
    %  Conversions between 3D Rotations. Modelling Simul. Mater. Sci. Eng.
    %  2015, 23 (8), 083501. https://doi.org/10.1088/0965-0393/23/8/083501.
    %--------------------------------------------------------------------------

    % qm = qmult(qBlab,qinv(qAlab),epsijk); %issue1: produces consistent results internally within GBdist(), but not during conversion, (passive convention)
    qm = qmult(qinv(qAlab), qBlab, epsijk); %issue2: produces consistent results in 5DOF-->octonion conversion, but not within GBdist(), (active convention)
end

function rnew = Lpr(p, r, epsijk)

    arguments
        p(:, 4) double {mustBeReal, mustBeFinite}
        r(:, 3) double {mustBeReal, mustBeFinite}
        epsijk(1, 1) double = 1
    end

    % LPR quaternion rotation active (epsijk==1) and passive (epsijk==-1)
    %--------------------------------------------------------------------------
    % Author(s): Sterling Baird
    % Date: 2020-07-27
    %
    % Inputs:
    %  qmis - misorientation quaternion
    %  r - boundary plane normal in cartesian coordinates
    %
    % Outputs:
    %  mA - rotated vector (i.e. r rotated by qmis)
    %
    % Usage:
    %  mA = Lpr(qmis,r);
    %
    % References:
    %  https://dx.doi.org/10.1088/0965-0393/23/8/083501 (eqn 24)
    %--------------------------------------------------------------------------
    %definitions
    p0 = p(:, 1);
    p = p(:, 2:4);
    pmag = vecnorm(p, 2, 2);

    %equation
    rnew = (p0.^2 - pmag.^2) .* r + 2 * dot(p, r, 2) .* p + 2 * epsijk * p0 .* cross(p, r, 2);

end
