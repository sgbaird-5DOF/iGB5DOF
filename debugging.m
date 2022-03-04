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

        if all(abs(projpts(:, ndim - ndegdim + 1:ndim)) < tol, 'all')
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
