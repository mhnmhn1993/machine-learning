% SLIC Simple Linear Iterative Clustering SuperPixels
%
% Implementation of Achanta, Shaji, Smith, Lucchi, Fua and Susstrunk's
% SLIC Superpixels
%
% Usage:   [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw)
%
% Arguments:  im - Image to be segmented.图像
%              k - Number of desired superpixels. Note that this is nominal
%                  the actual number of superpixels generated will generally
%                  be a bit larger, espiecially if parameter m is small.
%              m - Weighting factor between colour and spatial
%                  differences. Values from about 5 to 40 are useful.  Use a
%                  large value to enforce superpixels with more regular and
%                  smoother shapes. Try a value of 10 to start with.
%       seRadius - Regions morphologically smaller than this are merged with
%                  adjacent regions. Try a value of 1 or 1.5.  Use 0 to
%                  disable.
%         colopt - String 'mean' or 'median' indicating how the cluster
%                  colour centre should be computed. Defaults to 'mean'
%             mw - Optional median filtering window size.可选的中值滤波窗口大小
%                  Image compression can result in noticeable artifacts in the a*b* components of the image.  
%                  图像压缩可以导致在图像a*b*成分内明显的 人工产品 
%                  Median filtering can reduce this. mw can be a single value in which case the same median filtering is applied to each L* a* and b* components. 
%                  中值滤波可以减少。mw可以是一个包围相同median的应用于每个L* a* b*成分的中值滤波的单值
%                  Alternatively it can be a 2-vector where mw(1) specifies the median filtering window to be applied to L* and mw(2) is the median filtering window to be applied to a* and b*.
%                  作为一种选择，它可以是一个mw(1)指定中值滤波窗口的应用于L*and mw(2)的2维矢量，是应用于a* and b*的median滤波器窗口
% Returns:     l - Labeled image of superpixels. Labels range from 1 to
% k.标记的图像处理。标签范围从1到k
%             Am - Adjacency matrix of segments.  Am(i, j) indicates whether
%                  segments labeled i and j are connected/adjacent
%             Sp - Superpixel attribute structure array with fields:
%                   L  - Mean L* value
%                   a  - Mean a* value
%                   b  - Mean b* value
%                   r  - Mean row value
%                   c  - Mean column value
%                   stdL  - Standard deviation of L* 
%                   stda  - Standard deviation of a* 
%                   stdb  - Standard deviation of b* 
%                   N - Number of pixels
%                   edges - List of edge numbers that bound each
%                           superpixel. This field is allocated, but not set,
%                           by SLIC. Use SPEDGES for this.
%              d - Distance image giving the distance each pixel is from its
%                  associated superpixel centre.
%
% It is suggested that use of this function is followed by SPDBSCAN to perform a
% DBSCAN clustering of superpixels.  This results in a simple and fast
% segmentation of an image.
%这个SLIC函数后跟SPDBSCAN函数实现简单快速的图像分割
%
% Minor variations from the original algorithm as defined in Achanta et al's paper:
% 微小的变化,从原来Achanta的论文中的算法定义:
%
% - SuperPixel centres are initialised on a hexagonal grid rather than a square
%   one. This results in a segmentation that will be nominally 6-connected
%   which hopefully facilitates any subsequent post-processing that seeks to
%   merge superpixels.
% - Initial cluster positions are not shifted to point of lowest gradient
%   within a 3x3 neighbourhood because this will be rendered irrelevant the
%   first time cluster centres are updated.
% 超像素中心初始化一个六角网格,而不是一个正方形。
% 他名义上的分割结果将6-connected希望促进任何随后的后处理,旨在superpixels合并。
% 集群初始位置不转移到最低的梯度内3 x3附近,因为这将是无关紧要的第一次呈现集群中心更新。
% Reference: R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and
% S. Susstrunk. "SLIC Superpixels Compared to State-of-the-Art Superpixel
% Methods"  PAMI. Vol 34 No 11.  November 2012. pp 2274-2281.
%
% See also: SPDBSCAN, MCLEANUPREGIONS, REGIONADJACENCY, DRAWREGIONBOUNDARIES, RGB2LAB

% Copyright (c) 2013 Peter Kovesi
% Centre for Exploration Targeting
% School of Earth and Environment
% The University of Western Australia
% peter.kovesi at uwa edu au
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% Feb  2013
% July 2013 Super pixel attributes returned as a structure array

% Note that most of the computation time is not in the clustering, but rather
% in the region cleanup process.


function [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw, nItr, eim, We)
    
    if ~exist('colopt','var') || isempty(colopt), colopt = 'mean'; end
    if ~exist('mw','var')     || isempty(mw),         mw = 0;      end
    if ~exist('nItr','var')   || isempty(nItr),     nItr = 10;     end
    
    if exist('eim', 'var'), USEDIST = 0; else, USEDIST = 1; end
        
    MEANCENTRE = 1;
    MEDIANCENTRE = 2;
    
    if strcmp(colopt, 'mean')
        centre = MEANCENTRE;
    elseif strcmp(colopt, 'median')
        centre = MEDIANCENTRE;        
    else
        error('Invalid colour centre computation option');
    end
      [rows, cols, chan] = size(im);
    if chan ~= 3
        error('Image must be colour');
    end
  
    
    % Convert image to L*a*b* colourspace.  This gives us a colourspace that is
    % nominally perceptually uniform. This allows us to use the euclidean
    % distance between colour coordinates to measure differences between
    % colours.  Note the image becomes double after conversion.  We may want to
    % go to signed shorts to save memory.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 转换成L*a*b色彩空间。
    %这给了我们一个名义上的感知颜色空间。
    %这允许我们使用颜色坐标之间的欧氏距离来测量颜色之间的差异。
    %注意图像转换后变成double型。我们可能想去变成signed shorts型来节省内存。
    im = rgb2lab(im); 

    % Apply median filtering to colour components if mw has been supplied
    % and/or non-zero
    %如果mw被提供不是0，中值滤波应用于颜色组件
    if mw
        if length(mw) == 1
            mw(2) = mw(1);  % Use same filtering for L and chrominance
        end
        for n = 1:3
            im(:,:,n) = medfilt2(im(:,:,n), [mw(1) mw(1)]);
        end
    end
    
    % Nominal spacing between grid elements assuming hexagonal grid
    %在要素网格和假设的六边网格名义上的分隔
    S = sqrt(rows*cols / (k * sqrt(3)/2));
    
    % Get nodes per row allowing a half column margin at one end that alternates
    % from row to row
    nodeCols = round(cols/S - 0.5);
    % Given an integer number of nodes per row recompute S
    % 给一个每行节点的整数重新计算S
    S = cols/(nodeCols + 0.5); 

    % Get number of rows of nodes allowing 0.5 row margin top and bottom
    nodeRows = round(rows/(sqrt(3)/2*S));
    vSpacing = rows/nodeRows;

    % Recompute k
    k = nodeRows * nodeCols;
    
    % Allocate memory and initialise clusters, labels and distances.
    C = zeros(6,k);          % Cluster centre data  1:3 is mean Lab value,
                             % 4:5 is row, col of centre, 6 is No of pixels
    l = -ones(rows, cols);   % Pixel labels.
    d = inf(rows, cols);     % Pixel distances from cluster centres.
    
    % Initialise clusters on a hexagonal grid
    kk = 1;
    r = vSpacing/2;
    
    for ri = 1:nodeRows
        % Following code alternates the starting column for each row of grid
        % points to obtain a hexagonal pattern. Note S and vSpacing are kept
        % as doubles to prevent errors accumulating across the grid.
        if mod(ri,2), c = S/2; else, c = S;  end
        
        for ci = 1:nodeCols
            cc = round(c); rr = round(r);
            C(1:5, kk) = [squeeze(im(rr,cc,:)); cc; rr];
            c = c+S;
            kk = kk+1;
        end
        
        r = r+vSpacing;
    end
    
    % Now perform the clustering.  10 iterations is suggested but I suspect n
    % could be as small as 2 or even 1
    %迭代，可设置参数
    S = round(S);  % We need S to be an integer from now on
    
    for n = 1:nItr
       for kk = 1:k  % for each cluster

           % Get subimage around cluster
           rmin = max(C(5,kk)-S, 1);   rmax = min(C(5,kk)+S, rows); 
           cmin = max(C(4,kk)-S, 1);   cmax = min(C(4,kk)+S, cols); 
           subim = im(rmin:rmax, cmin:cmax, :);  
           assert(numel(subim) > 0)
           
           % Compute distances D between C(:,kk) and subimage
           if USEDIST
               D = dist(C(:, kk), subim, rmin, cmin, S, m);
           else
               D = dist2(C(:, kk), subim, rmin, cmin, S, m, eim, We);
           end

           % If any pixel distance from the cluster centre is less than its
           % previous value update its distance and label
           subd =  d(rmin:rmax, cmin:cmax);
           subl =  l(rmin:rmax, cmin:cmax);
           updateMask = D < subd;
           subd(updateMask) = D(updateMask);
           subl(updateMask) = kk;
           
           d(rmin:rmax, cmin:cmax) = subd;
           l(rmin:rmax, cmin:cmax) = subl;           
       end
       
       % Update cluster centres with mean values
       % 更新聚类中心与平均值
       C(:) = 0;
       for r = 1:rows
           for c = 1:cols
              tmp = [im(r,c,1); im(r,c,2); im(r,c,3); c; r; 1];
              C(:, l(r,c)) = C(:, l(r,c)) + tmp;
           end
       end
       
       % Divide by number of pixels in each superpixel to get mean values
       %在每个超像素除以数量的像素平均值
       for kk = 1:k 
           C(1:5,kk) = round(C(1:5,kk)/C(6,kk)); 
       end
       
       % Note the residual error, E, is not calculated because we are using a
       % fixed number of iterations 
       % 注意残留误差，E,不计算因为我们使用的是固定的迭代次数
    end
    
    % Cleanup small orphaned regions and 'spurs' on each region using
    % morphological opening on each labeled region.  The cleaned up regions are
    % assigned to the nearest cluster. The regions are renumbered and the
    % adjacency matrix regenerated.  This is needed because the cleanup is
    % likely to change the number of labeled regions.
    %清除小的孤立的区域，
    %重新编号和再生邻接矩阵这些区域。
    %这些需要被清理因为可能会改变标记区域的数量。
    if seRadius
        [l, Am] = mcleanupregions(l, seRadius);
    else
        l = makeregionsdistinct(l);
        [l, minLabel, maxLabel] = renumberregions(l);
        Am = regionadjacency(l);    
    end

    % Recompute the final superpixel attributes and write information into
    % the Sp struct array.
    %再计算最后超像素属性和信息写入Sp结构体数组。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = length(Am);
    Sp = struct('L', cell(1,N), 'a', cell(1,N), 'b', cell(1,N), ...
                'stdL', cell(1,N), 'stda', cell(1,N), 'stdb', cell(1,N), ...
                'r', cell(1,N), 'c', cell(1,N), 'N', cell(1,N));
    [X,Y] = meshgrid(1:cols, 1:rows);
    L = im(:,:,1);    
    A = im(:,:,2);    
    B = im(:,:,3);    
    for n = 1:N
        mask = l==n;
        nm = sum(mask(:));
        if centre == MEANCENTRE     
            Sp(n).L = sum(L(mask))/nm;
            Sp(n).a = sum(A(mask))/nm;
            Sp(n).b = sum(B(mask))/nm;
            
        elseif centre == MEDIANCENTRE
            Sp(n).L = median(L(mask));
            Sp(n).a = median(A(mask));
            Sp(n).b = median(B(mask));
        end
        
        Sp(n).r = sum(Y(mask))/nm;
        Sp(n).c = sum(X(mask))/nm;
        
        % Compute standard deviations of the colour components of each super
        % pixel. This can be used by code seeking to merge superpixels into
        % image segments.  Note these are calculated relative to the mean colour
        % component irrespective of the centre being calculated from the mean or
        % median colour component values.
        %计算每个超像素的颜色组件的标准差。
        %这可以通过代码去想办法合并超像素变为图像分割。
        %注意这些都是计算相对于中心的平均颜色组件无论从均值或中位数计算颜色组件值。
        Sp(n).stdL = std(L(mask));
        Sp(n).stda = std(A(mask));
        Sp(n).stdb = std(B(mask));

        Sp(n).N = nm;  % Record number of pixels in superpixel too.记录超像素中像素个数
    end
    
%-- dist -------------------------------------------
%
% Usage:  D = dist(C, im, r1, c1, S, m)
% 
% Arguments:   C - Cluster being considered 被考虑的簇
%             im - sub-image surrounding cluster centre 簇中心周围的子图
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing 网格步长
%              m - weighting factor between colour and spatial differences.权重
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre子图每个像素距离簇中心的距离
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance LAB颜色
% ds = sqrt(dx^2 + dy^2)         % Spatial distance 空间位置
%
% m is a weighting factor representing the nominal maximum colour distance
% expected so that one can rank colour similarity relative to distance
% similarity.  try m in the range [1-40] for L*a*b* space m为权重，范围为1-40
%
% ?? Might be worth trying the Geometric Mean instead ??
%  Distance = sqrt(dc * ds)
% but having a factor 'm' to play with is probably handy

% This code could be more efficient

function D = dist(C, im, r1, c1, S, m)

    % Squared spatial distance
    %    ds is a fixed 'image' we should be able to exploit this
    %    and use a fixed meshgrid for much of the time somehow...
    [rows, cols, chan] = size(im);
    [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
    x = x-C(4);  % x and y dist from cluster centre
    y = y-C(5);
    ds2 = x.^2 + y.^2;
    
    % Squared colour difference
    for n = 1:3
        im(:,:,n) = (im(:,:,n)-C(n)).^2;
    end
    dc2 = sum(im,3);
    
    D = sqrt(dc2 + ds2/S^2*m^2);
    
    
    
%--- dist2 ------------------------------------------
%
% Usage:  D = dist2(C, im, r1, c1, S, m, eim)
% 
% Arguments:   C - Cluster being considered
%             im - sub-image surrounding cluster centre
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing
%              m - weighting factor between colour and spatial differences.
%            eim - Edge strength sub-image corresponding to im 边缘强度
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% ds = sqrt(dx^2 + dy^2)         % Spatial distance
%
% m is a weighting factor representing the nominal maximum colour distance
% expected so that one can rank colour similarity relative to distance
% similarity.  try m in the range [1-40] for L*a*b* space
%

function D = dist2(C, im, r1, c1, S, m, eim, We)

    % Squared spatial distance
    %    ds is a fixed 'image' we should be able to exploit this
    %    and use a fixed meshgrid for much of the time somehow...
    [rows, cols, chan] = size(im);
    [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
    x = x-C(4);
    y = y-C(5);
    ds2 = x.^2 + y.^2;
    
    % Squared colour difference
    for n = 1:3
        im(:,:,n) = (im(:,:,n)-C(n)).^2;
    end
    dc2 = sum(im,3);
    
    % Combine colour and spatial distance measure
    D = sqrt(dc2 + ds2/S^2*m^2);
    
    % for every pixel in the subimage call improfile to the cluster centre
    % and use the largest value as the 'edge distance'
    rCentre = C(5)-r1;   % Cluster centre coords relative to this sub-image
    cCentre = C(4)-c1;
    de = zeros(rows,cols);
    for r = 1:rows
        for c = 1:cols
            v = improfile(eim,[c cCentre], [r rCentre]);
            de(r,c) = max(v);
        end
    end

    % Combine edge distance with weight, We with total Distance.边缘
    D = D + We * de;
    
