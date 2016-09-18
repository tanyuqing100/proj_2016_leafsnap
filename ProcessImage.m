%%% Leaf shape features extraction algorithms
%%% University of Saskatchewan, CMPT400
%%%
%%% Beatrice
%%% 2016

%%% File names and folders
datasize = 50; % for round-scanning model
samplesize = 1; % how many images to analyze from each species
dirpath_in = '160401in/'; % input images, 160401in/<species_folders>
dirpath_outimg = '160401outimg/'; % output images with a lot of markings
dirpath_outtxt = '160401outtxt/'; % output text data of concave-triangle model
dirpath_outplot = '160401outplot/'; % plot of text data

%species = dir(strcat(dirpath_in, 'metasequoia_glyptostroboide*'));
species = dir(strcat(dirpath_in, '*_*'));

%%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%
%%% Three output files for weka to use (done on 160209, 160312)
delete('weka1.txt');
output_weka1 = fopen('weka1.txt', 'a'); % round-scanning model
delete('weka2.txt');
output_weka2 = fopen('weka2.txt', 'a'); % 3-meta-properties model
delete('weka3.txt');
output_weka3 = fopen('weka3.txt', 'a'); % concave-triangle model

for output_weka = [output_weka1, output_weka2, output_weka3]
    fprintf(output_weka, '@relation leafdata\n');
    fprintf(output_weka, '@attribute ');
    fprintf(output_weka, '''imgname''');
    fprintf(output_weka, ' {');
    for sp = 1:length(species)
        %disp(species(sp).name);
        if (sp > 1) fprintf(output_weka, ','); end  
        for ss = 1:samplesize
            if (ss > 1) fprintf(output_weka, ','); end  
            fprintf(output_weka, '''%s%d''', species(sp).name, ss);
        end
    end
    fprintf(output_weka, '}\n');
end

for fn = 1:(datasize)
    fprintf(output_weka1, '@attribute ');
    fprintf(output_weka1, '''length%d''', fn);
    fprintf(output_weka1, ' numeric\n');
    fprintf(output_weka1, '@attribute ');
    fprintf(output_weka1, '''angle%d''', fn);
    fprintf(output_weka1, ' numeric\n');
end
fprintf(output_weka1, '@data\n');

fprintf(output_weka2, '@attribute '); fprintf(output_weka2, '''oblong'''); fprintf(output_weka2, ' numeric\n');
fprintf(output_weka2, '@attribute '); fprintf(output_weka2, '''fluffy'''); fprintf(output_weka2, ' numeric\n');
fprintf(output_weka2, '@attribute '); fprintf(output_weka2, '''jagged'''); fprintf(output_weka2, ' numeric\n');
fprintf(output_weka2, '@data\n');

%combine 2&3
%fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''oblong'''); fprintf(output_weka3, ' numeric\n');
%fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''fluffy'''); fprintf(output_weka3, ' numeric\n');
%fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''jagged'''); fprintf(output_weka3, ' numeric\n');
fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''numJagMain'''); fprintf(output_weka3, ' numeric\n');
fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''numJagTrivial'''); fprintf(output_weka3, ' numeric\n');
fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''peakRadiusAve'''); fprintf(output_weka3, ' numeric\n');
fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''peakRadiusStd'''); fprintf(output_weka3, ' numeric\n');
fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''peakBottomAve'''); fprintf(output_weka3, ' numeric\n');
fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''peakBottomStd'''); fprintf(output_weka3, ' numeric\n');
fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''curveBottonAve'''); fprintf(output_weka3, ' numeric\n');
fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''curveBottonStd'''); fprintf(output_weka3, ' numeric\n');
fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''topAngleAve'''); fprintf(output_weka3, ' numeric\n');
fprintf(output_weka3, '@attribute '); fprintf(output_weka3, '''sideAngleAve'''); fprintf(output_weka3, ' numeric\n');
fprintf(output_weka3, '@data\n');

%%% looping on species folders
for sp = 1:length(species)

ss = 1;
species_dir = strcat(dirpath_in, species(sp).name);
species_dir = strcat(species_dir, '/');
disp('species_dir');
disp(species_dir);
filename = dir(strcat(species_dir, '*.png'));
disp('filename');
disp(filename);

%%% looping on filenames in each species folder
for fn = 1:length(filename)

ss_str = sprintf('%d', ss);
input = strcat(species_dir,filename(fn).name);
BW = imread(input);
[B,L,N,A] = bwboundaries(BW);
figure; 
imshow(BW); 
hold on;

outfile = strcat(species(sp).name, ss_str);
outfile_txt = strcat(dirpath_outtxt, outfile);
outfile_txt = strcat(outfile_txt, '.txt');
delete(outfile_txt);
output_txt = fopen(outfile_txt, 'a');

fprintf(output_weka1, '''%s''', outfile);
fprintf(output_weka1, ', ');
fprintf(output_weka2, '''%s''', outfile);
fprintf(output_weka2, ', ');
fprintf(output_weka3, '''%s''', outfile);
fprintf(output_weka3, ', ');

%%% Adding text to images
text_lineSpacing = 20;
text_coordX = 10;
text_coordY = 20;
text_color = [0, 1, 1];

%%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%
%%% Getting info from regionprops() and boundary array
rp = regionprops(BW, 'MajorAxisLength', 'MinorAxisLength', ...
    'Area', 'FilledArea', 'ConvexArea', 'ConvexHull', ...
    'Perimeter', 'Centroid', 'Orientation', ...
    'Extrema', 'EulerNumber');
% all attributes are arrays since there are more than 1 objects in the image
rpArea_arr = cat(1, rp.FilledArea);
rpArea_sum = sum(rpArea_arr);
rpPerim_arr = cat(1,rp.Perimeter);
rpCentre_arr = cat(1, rp.Centroid);
rpOrient_arr = cat(1, rp.Orientation);
rpMajorlen_arr = cat(1, rp.MajorAxisLength);
rpMinorlen_arr = cat(1, rp.MinorAxisLength);
rpConvexhull_arr = cat(1,rp.ConvexHull);
rpConvexarea_arr = cat(1,rp.ConvexArea);

area_px2 = 0;
perim_sum = 0;
centreX_sum = 0;
centreY_sum = 0;
orient_sum = 0;
majorlen_sum = 0;
minorlen_sum = 0;
convexarea_px2 = 0;
boundary = [];

for obj = 1:length(rpArea_arr)
    if (rpArea_arr(obj,:) > rpArea_sum*0.01)
        % if an object occupies >1% of total area, it will not be ignored
        % each attribute is multiplied by area (weight)
        area_px2 = area_px2 + rpArea_arr(obj,:);
        perim_sum = perim_sum + rpPerim_arr(obj,:)*rpArea_arr(obj,:);
        centreX_sum = centreX_sum + rpCentre_arr(obj,1)*rpArea_arr(obj,:);
        centreY_sum = centreY_sum + rpCentre_arr(obj,2)*rpArea_arr(obj,:);
        orient_sum = orient_sum + rpOrient_arr(obj,:)*rpArea_arr(obj,:);
        majorlen_sum = majorlen_sum + rpMajorlen_arr(obj,:)*rpArea_arr(obj,:);
        minorlen_sum = minorlen_sum + rpMinorlen_arr(obj,:)*rpArea_arr(obj,:);
        convexarea_px2 = convexarea_px2 + rpConvexarea_arr(obj,:);
        disp(size(B{obj}'));
        disp(length(B{obj}'));
        boundary = [boundary, B{obj}'];
    end
end
boundary = boundary';

%%% Calculation on centroid, perimeter, orientation, axis length
radius_px = floor((area_px2/pi)^0.5);
perimeter_px = floor(perim_sum / area_px2);
centreX_px = floor(centreX_sum / area_px2);
centreY_px = floor(centreY_sum / area_px2);
orient_dgr = orient_sum / area_px2;
majorlen_px = majorlen_sum / area_px2;
minorlen_px = minorlen_sum / area_px2;

% print text
str = sprintf('area (px^2) = %d',area_px2); text(text_coordX,text_coordY,str,'Color',text_color); text_coordY = text_coordY + text_lineSpacing;
str = sprintf('perimeter (px) = %d',perimeter_px); text(text_coordX,text_coordY,str,'Color',text_color); text_coordY = text_coordY + text_lineSpacing;
str = sprintf('boundary size (px) = %d',length(boundary)); text(text_coordX,text_coordY,str,'Color',text_color); text_coordY = text_coordY + text_lineSpacing;
str = sprintf('average radius (px) = %d',radius_px); text(text_coordX,text_coordY,str,'Color',text_color); text_coordY = text_coordY + text_lineSpacing;
str = sprintf('centroid (px,px) = (%d,%d)',centreX_px,centreY_px); text(text_coordX,text_coordY,str,'Color',text_color); text_coordY = text_coordY + text_lineSpacing;

% draw out centroid and axis in red
plot(centreX_px, centreY_px, 'r*');
axisX = centreX_px - majorlen_px/2 * cosd(orient_dgr); axisY = centreY_px - majorlen_px/2 * sind(orient_dgr); plot(axisX, axisY, 'r*');
axisX = centreX_px + majorlen_px/2 * cosd(orient_dgr); axisY = centreY_px + majorlen_px/2 * sind(orient_dgr); plot(axisX, axisY, 'r*');
axisX = centreX_px - minorlen_px/2 * cosd(orient_dgr+90); axisY = centreY_px - minorlen_px/2 * sind(orient_dgr+90); plot(axisX, axisY, 'r*');
axisX = centreX_px + minorlen_px/2 * cosd(orient_dgr+90); axisY = centreY_px + minorlen_px/2 * sind(orient_dgr+90); plot(axisX, axisY, 'r*');
for ml = 1:majorlen_px
    axisX = centreX_px + (ml - majorlen_px/2) * cosd(orient_dgr);
    axisY = centreY_px + (ml - majorlen_px/2) * sind(orient_dgr);
    plot(axisX, axisY, '-r');
end
for ml = 1:minorlen_px
    axisX = centreX_px + (ml - minorlen_px/2) * cosd(orient_dgr+90);
    axisY = centreY_px + (ml - minorlen_px/2) * sind(orient_dgr+90);
    plot(axisX, axisY, '-r');
end

%%% 3-meta-properties model features (done on 160209)
% overall shape
oblong = minorlen_px / majorlen_px; % the smaller the more oblong
jagged = (area_px2*pi*4) / (perimeter_px^2); % the smaller the more jagged
fluffy = area_px2 / (majorlen_px*minorlen_px); % the smaller the more fluffy
%fluffy = area_px2 / convexarea_px2;
fprintf(output_weka2, '%f, %f, %f\n', oblong, fluffy, jagged);

%combine 1&2
%fprintf(output_weka3, '%f, %f, %f, ', oblong, fluffy, jagged);

%%% concave-triangle model features (done on 160308)
% arithmetic on convex hull
% (a) number of jags; (b) peakHeight/radius (ave and std);
% (c) peakHeight/bottom; (d) curveLength/bottom; (e) solve triangle;
jagMain = []; % get a number of features from main jags
numJagMain = 0;
numJagTrivial = 0; % only count the number
fprintf('Matching convex points with boundary index...\n');
convexIndex_arr = [];
for bd=1:length(boundary)
    for ch=1:length(rpConvexhull_arr)
        if ( (boundary(bd,2) < rpConvexhull_arr(ch,1)+1 && boundary(bd,2) > rpConvexhull_arr(ch,1)-1) ...
            && (boundary(bd,1) < rpConvexhull_arr(ch,2)+1 && boundary(bd,1) > rpConvexhull_arr(ch,2)-1) )
            fprintf('index=%d bdX=%d bdY=%d chX=%f chy=%f\n', ...
                bd, boundary(bd,2), boundary(bd,1), rpConvexhull_arr(ch,1), rpConvexhull_arr(ch,2));
            convexIndex_arr = [convexIndex_arr, bd];
            plot(boundary(bd,2), boundary(bd,1), 'c*');
            %for i=1:10
            %    BW(boundary(bd,1)+i, boundary(bd,2)+i) = 0;
            %end
            %fprintf('%f %f\n', BW(boundary(bd,2), boundary(bd,1)), BW(boundary(bd,2)+2, boundary(bd,1)+2));
            break;
        end
    end
end
fprintf('\n');

fprintf('Solving concave triangles...\n');
convexX1 = 0;
convexY1 = 0;
convexX2 = 0;
convexY2 = 0;
ch_prev = 0;
triangData = [];
triangThres1 = floor(radius_px/10); % interpolate point
triangThres2 = floor(radius_px/20); % height smaller than this will be considered flat area
for ch=1:length(convexIndex_arr)
    convexX1 = boundary(convexIndex_arr(ch),2);
    convexY1 = boundary(convexIndex_arr(ch),1);
    if ((convexX1 ~= 0) && (convexY1 ~= 0) && (convexX2 ~= 0) && (convexY2 ~= 0))
        fprintf('index=%d convexX1=%d convexY1=%d convexX2=%d convexY2=%d\n', convexIndex_arr(ch), convexX1, convexY1, convexX2, convexY2);
        triangBottom = ((convexX1 - convexX2)^2 + (convexY1 - convexY2)^2)^0.5;
        % line equation: paramA * x + param * y + 1000 = 0
        paramA = (convexY1 - convexY2) / (convexX1 - convexX2);
        paramB = -1;
        paramC = convexY1 - paramA * convexX1;
        paramA = paramA / paramC * 1000;
        paramB = paramB / paramC * 1000;
        paramC = paramC / paramC * 1000;
        fprintf('paramA=%f paramB=%f paramC=%d\n', paramA, paramB, paramC);
        if isnan(paramA) || isnan(paramB) || isnan(paramC)
            continue;
        end
        
        centerSide = paramA * centreX_px + paramB * centreY_px + paramC;
        localAbsolute = []; %(minmax, index, height, cumulatedCurve)
        numAbsolute = 0;
        disp(size(localAbsolute));
        localMax = [0 0.0]; %(index, height)
        triangArea = 0;
        curveDist = 0;
        convexX4 = 0;
        convexY4 = 0;
        for point=convexIndex_arr(ch_prev):convexIndex_arr(ch)
            convexX3 = boundary(point,2);
            convexY3 = boundary(point,1);
            
            % increment curve length
            if ((convexX3 ~= 0) && (convexY3 ~= 0) && (convexX4 ~= 0) && (convexY4 ~= 0))
                stepDist = ((convexX3 - convexX4)^2 + (convexY3 - convexY4)^2)^0.5;
                curveDist = curveDist + stepDist;
            end
            
            % relation to centroid
            pointHeight = (paramA * convexX3 + paramB * convexY3 + paramC) / (paramA^2 + paramB^2)^0.5;
            if ( (pointHeight * centerSide > 0) )
                % concave; same side with centroid
                pointHeight = abs(pointHeight);
            else
                % convex; opposite side with centroid
                pointHeight = -abs(pointHeight);
            end
            triangArea = triangArea + pointHeight;
            
            % local maximum or minimum
            interpIndex_arr = zeros(1);
            interpDist_arr = zeros(1);
            for i=1:5
                interpIndex_arr(i) = point-(6-i);
                if interpIndex_arr(i) < 1
                    interpIndex_arr(i) = interpIndex_arr(i)+length(boundary);
                end
            end
            for i=6:10
                interpIndex_arr(i) = point+(i-5);
                if interpIndex_arr(i) > length(boundary)
                    interpIndex_arr(i) = interpIndex_arr(i)-length(boundary);
                end
            end
            peakCheck = 1;
            troughCheck = 1;
            special = 0;
            for i=1:10
                interpDist_arr(i) = (paramA * boundary(interpIndex_arr(i),2) + ...
                    paramB * boundary(interpIndex_arr(i),1) + paramC) / (paramA^2 + paramB^2)^0.5;
                if ( (interpDist_arr(i) * centerSide > 0) )
                    interpDist_arr(i) = abs(interpDist_arr(i));
                else
                    interpDist_arr(i) = -abs(interpDist_arr(i));
                end
                if interpDist_arr(i) > pointHeight
                    peakCheck = peakCheck * 0;
                    troughCheck = troughCheck * 1;
                end
                if interpDist_arr(i) < pointHeight
                    peakCheck = peakCheck * 1;
                    troughCheck = troughCheck * 0;
                end
            end
            if peakCheck == 1
                fprintf('peak at %d, %f\n', point, pointHeight);
                if pointHeight > localMax(:,2)
                    localMax(:,1) = point;
                    localMax(:,2) = pointHeight;
                end
                numAbsolute = numAbsolute + 1;
                localAbsolute(numAbsolute,:) = [1 point pointHeight curveDist];
                special = 1;
            end
            if troughCheck == 1
                fprintf('trough at %d, %f\n', point, abs(pointHeight));
                numAbsolute = numAbsolute + 1;
                localAbsolute(numAbsolute,:) = [-1 point pointHeight curveDist];
                %special = 1;
            end
            
            fprintf(output_txt, '%d %f %d\n', point, pointHeight, special);
            convexX4 = convexX3;
            convexY4 = convexY3;
        end
        
        aveHeight = triangArea / triangBottom;
        fprintf('aveHeight = %f\n', aveHeight );
        fprintf('localMax = %d %f\n\n', localMax(:,1), localMax(:,2));
        
        if aveHeight < 3
            convexX2 = convexX1;
            convexY2 = convexY1;
            ch_prev = ch;
            fprintf('flat area, skip checking\n');
            continue;
        end
        
        disp(localAbsolute);
        
        % situation where more than 1 maximums
        maxIndex = 1;
        minIndex = 1;
        for ab=1:numAbsolute
            if localAbsolute(ab,2) == localMax(:,1)
                maxIndex = ab;
                localAbsolute(ab,1) = 2;
                break;
            end
        end
        search = -1; % -1 for searching troughs
        ab=maxIndex;
        while ab>=1
            %fprintf('maxIndex=%d minIndex=%d ab=%d diff=%f\n', maxIndex,minIndex, ab, localAbsolute(maxIndex,3)-localAbsolute(ab,3));
            if search == -1 && localAbsolute(ab,1) == -1 && ab ~= 1 && ...
            localAbsolute(maxIndex,3)-localAbsolute(ab,3) > aveHeight*0.5 
                localAbsolute(ab,1) = -2;
                minIndex = ab;
                search = 1; % 1 for searching peaks
            end
            if search == 1 && localAbsolute(ab,1) == 1 && ...
            localAbsolute(ab,3)-localAbsolute(minIndex,3)> aveHeight*0.5 
                localAbsolute(ab,1) = 2;
                minIndex = 0;
                search = -1;
            end
            ab = ab-1;
        end
        if minIndex ~= 0
            localAbsolute(minIndex,1) = -1;
            minIndex = 0;
        end
        search = -1; % -1 for searching troughs
        ab=maxIndex;
        disp(length(localAbsolute));
        disp(size(localAbsolute));
        disp(numAbsolute);
        while ab<=numAbsolute
            %fprintf('maxIndex=%d minIndex=%d ab=%d diff=%f\n', maxIndex,minIndex, ab, localAbsolute(maxIndex,3)-localAbsolute(ab,3));
            if search == -1 && localAbsolute(ab,1) == -1 && ab ~= numAbsolute && ...
            localAbsolute(maxIndex,3)-localAbsolute(ab,3) > aveHeight*0.5 
                localAbsolute(ab,1) = -2;
                minIndex = ab;
                search = 1; % 1 for searching peaks
            end
            if search == 1 && localAbsolute(ab,1) == 1 && ...
            localAbsolute(ab,3)-localAbsolute(minIndex,3) > aveHeight*0.5 
                localAbsolute(ab,1) = 2;
                minIndex = 0;
                search = -1;
            end
            ab = ab+1;
        end
        if minIndex ~= 0
            localAbsolute(minIndex,1) = -1;
            minIndex = 0;
        end
        disp(localAbsolute);
        
        triangStart = convexIndex_arr(ch_prev);
        triangEnd = 0;
        curveStart = 0;
        curveEnd = 0;
        maxIndex = 0;
        peakIndex = 0;
        for ab=1:numAbsolute
            if localAbsolute(ab,1) == 1
                numJagTrivial=numJagTrivial + 1;
            end
            if localAbsolute(ab,1) == 2
                maxIndex = ab;
                peakIndex = localAbsolute(maxIndex,2);
                peakHeight = localAbsolute(maxIndex,3);
            end
            if localAbsolute(ab,1) == -2
                triangEnd = localAbsolute(ab,2);
                curveEnd = localAbsolute(ab,4);
            end
            if ab == numAbsolute
                triangEnd = convexIndex_arr(ch);
                curveEnd = curveDist;
            end
            
            if triangEnd ~= 0 && peakIndex ~= 0
                triangBottom = ((boundary(triangStart,2) - boundary(triangEnd,2))^2 + (boundary(triangStart,1) - boundary(triangEnd,1))^2)^0.5;
                peakLeftDist = ((boundary(peakIndex,2) - boundary(triangStart,2))^2 + (boundary(peakIndex,1) - boundary(triangStart,1))^2)^0.5;
                peakRightDist = ((boundary(peakIndex,2) - boundary(triangEnd,2))^2 + (boundary(peakIndex,1) - boundary(triangEnd,1))^2)^0.5;
                circumAngle = acosd( (peakLeftDist^2 + peakRightDist^2 - triangBottom^2)/(2*peakLeftDist*peakRightDist) );
                
                leftMidPoint = floor((peakIndex + triangStart)/2);
                rightMidPoint = floor((peakIndex + triangEnd)/2);
                leftMidPointLD = ((boundary(leftMidPoint,2) - boundary(triangStart,2))^2 + (boundary(leftMidPoint,1) - boundary(triangStart,1))^2)^0.5;
                leftMidPointRD = ((boundary(leftMidPoint,2) - boundary(triangEnd,2))^2 + (boundary(leftMidPoint,1) - boundary(triangEnd,1))^2)^0.5;
                rightMidPointLD = ((boundary(rightMidPoint,2) - boundary(triangStart,2))^2 + (boundary(rightMidPoint,1) - boundary(triangStart,1))^2)^0.5;
                rightMidPointRD = ((boundary(rightMidPoint,2) - boundary(triangEnd,2))^2 + (boundary(rightMidPoint,1) - boundary(triangEnd,1))^2)^0.5;
                centralAngleL = acosd( (leftMidPointLD^2 + leftMidPointRD^2 - triangBottom^2)/(2*leftMidPointLD*leftMidPointRD) );
                centralAngleR = acosd( (rightMidPointLD^2 + rightMidPointRD^2 - triangBottom^2)/(2*rightMidPointLD*rightMidPointRD) );
                
                numJagMain=numJagMain+1;
                jagMain(numJagMain,:) = [...
                    triangStart peakIndex triangEnd ...
                    peakHeight/radius_px ...
                    peakHeight/triangBottom ...
                    (curveEnd-curveStart)/triangBottom ...
                    circumAngle ...
                    (centralAngleL+centralAngleR)/2 ...
                    ];
                triangStart = triangEnd;
                triangEnd = 0;
                curveStart = curveEnd;
                peakIndex = 0;
            end
        end
    end
    convexX2 = convexX1;
    convexY2 = convexY1;
    ch_prev = ch;
end
fprintf('\n');
fprintf('numJagTrivial=%d\n', numJagTrivial);
disp('jagMain');
disp(jagMain);

if isempty(jagMain)
    fprintf(output_weka3, '0, 0, 0, 0, 0, 0, 0, 0, 0, 0\n');
    %fprintf(output_weka3, '0, 0, 0, 0, 0, 0, 0, 0\n');
else
    peakRadiusAve = mean(jagMain(:,4));
    peakRadiusStd = std(jagMain(:,4));
    peakBottomAve = mean(jagMain(:,5));
    peakBottomStd = std(jagMain(:,5));
    curveBottonAve = mean(jagMain(:,6));
    curveBottonStd = std(jagMain(:,6));
    circumAngleAve = mean(jagMain(:,7));
    centralAngleAve = mean(jagMain(:,8));

    fprintf('%d, %d, %f, %f, %f, %f, %f, %f, %f, %f\n', numJagMain, numJagTrivial, ...
        peakRadiusAve, peakRadiusStd, peakBottomAve, peakBottomStd, ...
        curveBottonAve, curveBottonStd, circumAngleAve, centralAngleAve);

    fprintf(output_weka3, '%d, %d, %f, %f, %f, %f, %f, %f, %f, %f\n', numJagMain, numJagTrivial, ...
        peakRadiusAve, peakRadiusStd, peakBottomAve, peakBottomStd, ...
        curveBottonAve, curveBottonStd, circumAngleAve, centralAngleAve);

    %{
    %combine 2&3
    fprintf(output_weka3, '%d, %d, %f, %f, %f, %f, %f, %f\n', numJagMain, numJagTrivial, ...
        peakRadiusAve, peakRadiusStd, peakBottomAve, peakBottomStd, ...
        curveBottonAve, curveBottonStd);
    %}
end

%%% start of round-scanning model features
%%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%
%%% Clip shape detection (marking out stem) (done on 1512__)
%%%     detecting clip shape by walking along boundary points and measuring
%%%     distance
%%%     time: O(n^2)

clipThres_const = 20;
clipStart_index = 0;
clipEnd_index = 0;
clipCount_num = 0;
clip_dist = 0;
test_dist = 0;

clipThres = floor(radius_px/clipThres_const);
clip_intv1 = floor(clipThres*2.5);
clip_intv2 = floor(clipThres*7.5);
str = sprintf('clip threshold (px) = %d',clipThres); text(text_coordX,text_coordY,str,'Color',text_color); text_coordY = text_coordY + text_lineSpacing;
bd = 1;
while (bd <= length(boundary))
    if (bd > clip_intv1)
        test_disty = boundary(bd,1) - boundary(bd-clip_intv1,1);
        test_distx = boundary(bd,2) - boundary(bd-clip_intv1,2);
    else
        test_disty = boundary(bd,1) - boundary(bd-clip_intv1+length(boundary),1);
        test_distx = boundary(bd,2) - boundary(bd-clip_intv1+length(boundary),2);
    end
    test_dist = (test_disty^2 + test_distx^2)^0.5;
    
    if (test_dist < clipThres)
        clip_dist = test_dist;
        p1 = 1;
        p2 = -1;
        
        % decide size of the clip shape, mark out in grey
        while (clip_dist < clipThres)
            q1 = bd;
            q2 = bd - clip_intv1;
            if (q1+p1 > length(boundary)) q1 = q1 - length(boundary); end
            if (q2+p2 < 1) q2 = q2 + length(boundary); end
            clip_disty1 = boundary(q1+p1,1) - boundary(q2+p2,1);
            clip_disty2 = boundary(q1+p1+1,1) - boundary(q2+p2,1);
            clip_disty3 = boundary(q1+p1,1) - boundary(q2+p2-1,1);
            clip_distx1 = boundary(q1+p1,2) - boundary(q2+p2,2);
            clip_distx2 = boundary(q1+p1+1,2) - boundary(q2+p2,2);
            clip_distx3 = boundary(q1+p1,2) - boundary(q2+p2-1,2);
            clip_dist1 = (clip_disty1^2 + clip_distx1^2)^0.5;
            clip_dist2 = (clip_disty2^2 + clip_distx2^2)^0.5;
            clip_dist3 = (clip_disty3^2 + clip_distx3^2)^0.5;
            if (clip_dist2 < clip_dist1)
               if (clip_dist3 < clip_dist2)
                   p2 = p2 - 1;
                   clip_dist = clip_dist3;
               else
                   p1 = p1 + 1;
                   clip_dist = clip_dist2;
               end
            else
               if (clip_dist3 < clip_dist1)
                   p2 = p2 - 1;
                   clip_dist = clip_dist3;
               else
                   clip_dist = clip_dist1;
               end
            end
            plot(boundary(q1+p1,2), boundary(q1+p1,1),'Color',[0.5,0.5,0.5],'Marker','.');
            plot(boundary(q2+p2,2), boundary(q2+p2,1),'Color',[0.5,0.5,0.5],'Marker','.');
            p1 = p1 + 1;
            p2 = p2 - 1;
        end
        
        p3 = p1 + clipThres;
        p4 = p2 - clipThres;
        if (q1+p3 > length(boundary)) p3 = p3 - length(boundary); end
        if (q2+p4 < 1) p4 = p4 + length(boundary); end
        clip_disty4 = boundary(q1+p3,1) - boundary(q2+p4,1);
        clip_distx4 = boundary(q1+p3,2) - boundary(q2+p4,2);
        clip_dist4 = (clip_disty4^2 + clip_distx4^2)^0.5;
        clip_disty5 = boundary(q1+p3,1) - boundary(bd,1);
        clip_distx5 = boundary(q1+p3,2) - boundary(bd,2);
        clip_dist5 = (clip_disty5^2 + clip_distx5^2)^0.5;
        clip_disty6 = boundary(q2+p4,1) - boundary(bd,1);
        clip_distx6 = boundary(q2+p4,2) - boundary(bd,2);
        clip_dist6 = (clip_disty6^2 + clip_distx6^2)^0.5;
        %plot(boundary(q1+p3,2), boundary(q1+p3,1),'y*');
        %plot(boundary(q2+p4,2), boundary(q2+p4,1),'y*');
        
        if (q1+p1 > q2+p2)
            clip_dist7 = (q1+p1) - (q2+p2);
        else
            clip_dist7 = (q1+p1) - (q2+p2) + length(boundary);
        end
        
        % decide if it's qualified to be a stem
        if (((clip_dist7 > clip_intv2) && (clip_dist4 > clip_intv1)) || ...
            ((clip_dist7 > clip_intv2*0.5) && (clip_dist5 + clip_dist6 < clip_dist7*0.7)) )
            if (q1+p1 > bd)
                %str = sprintf('q1+p1=%f, j=%f', q1+p1, j);
                %text(text_coordX,text_coordY,str,'Color',text_color);
                %text_coordY = text_coordY + text_lineSpacing;
                bd = q1+p1;
            else
                %str = sprintf('q1+p1=%f, j=%f', q1+p1, j);
                %text(text_coordX,text_coordY,str,'Color',text_color);
                %text_coordY = text_coordY + text_lineSpacing;
                bd = length(boundary);
            end
            if (clipCount_num == 0)
                clipEnd_index = q1+p1;
                clipStart_index = q2+p2;
                clipCount_num = clipCount_num + 1;
            elseif (clipCount_num > 0)
                % if more than one stem, use the one closest to centroid
                % (botanically easier for water transmission)
                if ((boundary(q1+p1,2)-centreX_px)^2+(boundary(q1+p1,1)-centreY_px)^2 < ...
                    (boundary(clipEnd_index,2)-centreX_px)^2+(boundary(clipEnd_index,1)-centreY_px)^2 )
                    clipEnd_index = q1+p1;
                    clipStart_index = q2+p2;
                end                
            end
        end
    end
    bd = bd + 1;
end
str = sprintf('number of stems = %d',clipCount_num);
text(text_coordX,text_coordY,str,'Color',text_color);
text_coordY = text_coordY + text_lineSpacing;

%%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%
%%% Symmetry detection (if failed to detect any stem) (done on 1512__)
%%%     detecting reflective symmetry by measuring angles and distances
%%%     between interpolating points
%%%     time: O(n^2)

symthres_param = 30;
lineangle_ratio = 0;
sym_index = 0;
stem_index = 0;
stem_degr = 0;

sym_intv = floor(length(boundary)/symthres_param);
sym_halfsize = 6;
offset_array = zeros(sym_halfsize*2);
angle_array = zeros(sym_halfsize*2);
lineseg_array = zeros(sym_halfsize);
if (clipCount_num == 0)
    for bd = 1:length(boundary)
        angle_sum = 0; % as small as possible
        line_sum = 0; % as large as possible
        for h = 1:sym_halfsize
            offset_array(h) = bd - sym_intv*(sym_halfsize+1-h);
            if (offset_array(h) < 1)
                offset_array(h) = offset_array(h) + length(boundary); end
            offset_array(sym_halfsize*2+1-h) = bd + sym_intv*(sym_halfsize+1-h);
            if (offset_array(sym_halfsize*2+1-h) > length(boundary)) 
                offset_array(sym_halfsize*2+1-h) = offset_array(sym_halfsize*2+1-h) - length(boundary); end
            lineseg_array(h) = ( ...
                (boundary(offset_array(h),1)-boundary(offset_array(sym_halfsize*2+1-h),1))^2 + ...
                (boundary(offset_array(h),2)-boundary(offset_array(sym_halfsize*2+1-h),2))^2 )^0.5;
            line_sum = line_sum + lineseg_array(h);
        end

        for h = 1:(sym_halfsize-1)
            angle_array(h) = atand( ...
                (boundary(offset_array(h+1),1)-boundary(offset_array(h),1)) / ...
                (boundary(offset_array(h+1),2)-boundary(offset_array(h),2)) );
            angle_array(sym_halfsize*2+1-h) = atand( ...
                (boundary(offset_array(sym_halfsize*2-h),1)-boundary(offset_array(sym_halfsize*2+1-h),1)) / ...
                (boundary(offset_array(sym_halfsize*2-h),2)-boundary(offset_array(sym_halfsize*2+1-h),2)) );                
        end
        angle_array(sym_halfsize) = atand( ...
            (boundary(bd,1)-boundary(offset_array(sym_halfsize),1)) / ...
            (boundary(bd,2)-boundary(offset_array(sym_halfsize),2)) );
        angle_array(sym_halfsize+1) = atand( ...
            (boundary(bd,1)-boundary(offset_array(sym_halfsize+1),1)) / ...
            (boundary(bd,2)-boundary(offset_array(sym_halfsize+1),2)) );                

        for h = 1:(sym_halfsize-1)
            angle_diff1 = angle_array(h+1) - angle_array(h);
            angle_diff2 = angle_array(sym_halfsize*2-h) - angle_array(sym_halfsize*2+1-h);
            angle_sum = angle_sum + (angle_diff1 + angle_diff2)^2;
        end
        
        % self invented formula for symmetry detection
        if ((line_sum^3/length(boundary)) / (angle_sum) > lineangle_ratio || lineangle_ratio == 0)
            sym_index = bd;
            lineangle_ratio = (line_sum^3/length(boundary)) / (angle_sum);
        end
    end
    clipStart_index = sym_index;
    clipEnd_index = sym_index;
end

if (clipStart_index > floor(length(boundary)/2))
    clipstart_oppo = clipStart_index - floor(length(boundary)/2);
else
    clipstart_oppo = clipStart_index + floor(length(boundary)/2);
end
if (clipEnd_index > floor(length(boundary)/2))
    clipend_oppo = clipEnd_index - floor(length(boundary)/2);
else
    clipend_oppo = clipEnd_index + floor(length(boundary)/2);
end
%clip_center_dist = ...
clip_center_dist1 = ((boundary(clipStart_index,2)-centreX_px)^2 + (boundary(clipStart_index,1)-centreY_px)^2)^0.5;
clip_center_dist2 = ((boundary(clipEnd_index,2)-centreX_px)^2 + (boundary(clipEnd_index,1)-centreY_px)^2)^0.5;
%clipreverse_center_dist = ...
clipreverse_center_dist1 = ((boundary(clipstart_oppo,2)-centreX_px)^2 + (boundary(clipstart_oppo,1)-centreY_px)^2)^0.5;
clipreverse_center_dist2 = ((boundary(clipend_oppo,2)-centreX_px)^2 + (boundary(clipend_oppo,1)-centreY_px)^2)^0.5;

if (clipStart_index ~= clipEnd_index || clip_center_dist1 < clipreverse_center_dist1 || clip_center_dist2 < clipreverse_center_dist2)
    plot(boundary(clipStart_index,2), boundary(clipStart_index,1),'b*','MarkerSize',10);
    plot(boundary(clipEnd_index,2), boundary(clipEnd_index,1),'b*','MarkerSize',10);
    stem_index = clipEnd_index;
else
    clip_center_slope = (boundary(clipStart_index,1)-centreY_px) / (boundary(clipStart_index,2)-centreX_px);
    clip_center_intersect = centreY_px - centreX_px*clip_center_slope;
    
    for n = 1:(length(boundary)/5)
        if (clipstart_oppo+n > length(boundary))
            if ((floor(boundary(clipstart_oppo+n-length(boundary),1)-boundary(clipstart_oppo+n-length(boundary),2)*clip_center_slope-clip_center_intersect))^2 <= 10)
                stem_index = clipstart_oppo+n;
                break;
            end
        else
            if ((floor(boundary(clipstart_oppo+n,1)-boundary(clipstart_oppo+n,2)*clip_center_slope-clip_center_intersect))^2 <= 10)
                stem_index = clipstart_oppo+n;
                break;
            end
        end
        if (clipstart_oppo-n < 1)
            if ((floor(boundary(clipstart_oppo-n+length(boundary),1)-boundary(clipstart_oppo-n+length(boundary),2)*clip_center_slope-clip_center_intersect))^2 <= 10)
                stem_index = clipstart_oppo-n;
                break;
            end
        else
            if ((floor(boundary(clipstart_oppo-n,1)-boundary(clipstart_oppo-n,2)*clip_center_slope-clip_center_intersect))^2 <= 10)
                stem_index = clipstart_oppo-n;
                break;
            end
        end
    end
    
    if (stem_index ~= 0)
        while (stem_index > length(boundary) || stem_index < 1)
            if (stem_index > length(boundary))
                stem_index = stem_index - length(boundary);
            end
            if (stem_index < 1)
                stem_index = stem_index + length(boundary);
            end
        end
        plot(boundary(stem_index,2), boundary(stem_index,1),'g*','MarkerSize',10);
    else
        stem_index = clipEnd_index;
    end
    plot(boundary(stem_index,2), boundary(stem_index,1),'b*','MarkerSize',10);
end

%%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%% %%%%%%%%%%%%
%%%(f) Collect data into txt (finally) while plotting out

stem_degr = atand((boundary(stem_index,1)-centreY_px) / (boundary(stem_index,2)-centreX_px));
str = sprintf('root position (index) = %d',stem_index);
text(text_coordX,text_coordY,str,'Color',text_color);
text_coordY = text_coordY + text_lineSpacing;
str = sprintf('root axis (deg) = %f',stem_degr);
text(text_coordX,text_coordY,str,'Color',text_color);
text_coordY = text_coordY + text_lineSpacing;

%centreY_px = boundary(stem_index,1);
%centreX_px = boundary(stem_index,2);
i = 0; % data line counter
j = 1; % now length(boundary) counter
d = length(boundary) / datasize;
if (clipStart_index > clipEnd_index)
    d = (clipStart_index - clipEnd_index) / datasize;
else
    d = (length(boundary) + clipStart_index - clipEnd_index) / datasize;
end
disp('d =');
disp(d);
ctr_disty = boundary(j,1) - centreY_px;
ctr_distx = boundary(j,2) - centreX_px;
ctr_dist = (ctr_disty^2 + ctr_distx^2)^0.5;
curangle_degr = atand(ctr_disty / ctr_distx);
anglediff_degr = curangle_degr - stem_degr;
anglediffprev_degr = anglediff_degr;

% draw out new identified axis in blue
for m = 1:ctr_dist
    axisx = centreX_px - m * cosd(curangle_degr);
    axisy = centreY_px - m * sind(curangle_degr);
    plot(axisx, axisy, '-b');
    axisx = centreX_px + m * cosd(curangle_degr);
    axisy = centreY_px + m * sind(curangle_degr);
    plot(axisx, axisy, '-b');
end

while (i < datasize)
    i = i + 1;
    j = j + d;
    while (floor(j) > length(boundary) || floor(j) < 1)
        if (floor(j) > length(boundary))
            j = j - length(boundary);
        end
        if (floor(j) < 1)
            j = j + length(boundary);
        end
    end
    if (i > 1) fprintf(output_weka1, ', '); end

    ctr_disty = boundary(floor(j),1) - centreY_px;
    ctr_distx = boundary(floor(j),2) - centreX_px;
    ctr_dist = (ctr_disty^2 + ctr_distx^2)^0.5;
    curangle_degr = atand(ctr_disty / ctr_distx);
    anglediff_degr = curangle_degr - stem_degr;
    if (anglediff_degr < 0)
        anglediff_degr = anglediff_degr + 180;
    end
    while (anglediff_degr < anglediffprev_degr - 90)
        %fprintf(output, 'anglediff_degr=%f anglediffprev_degr=%f\n', anglediff_degr, anglediffprev_degr);
        anglediff_degr = anglediff_degr + 180;
    end
    anglediffprev_degr = anglediff_degr;

    %data print 1
    %fprintf(output_txt, '%d %d %f %f\n', i, floor(j), ctr_dist/radius_px, anglediff_degr);
    
    %data print 2
    fprintf(output_weka1, '%f, %f', ctr_dist/radius_px, anglediff_degr);
    
    % quadrants are different from drawing (y is upsidedown), 
    % marking with colors
    xx = boundary(floor(j),2);
    yy = boundary(floor(j),1);
    if (ctr_disty > 0)
        if (ctr_distx > 0)
            plot(xx, yy,'y.'); %r
            %fprintf(output_txt, 'r\n');
        else
            plot(xx, yy,'y.'); %y
            %fprintf(output_txt, 'y\n');
        end
    else
        if (ctr_distx < 0)
            plot(xx, yy,'y.'); %c
            %fprintf(output_txt, 'c\n');
        else
            plot(xx, yy,'y.'); %b
            %fprintf(output_txt, 'b\n');
        end
    end
end

% save the images with a lot of markings on it, into /outimg
fprintf(output_weka1, '\n');
outfile_img = strcat(dirpath_outimg,outfile);
outfile_img = strcat(outfile_img,'.png');
print(outfile_img,'-dpng')

% read text files in /outtxt and plot images into /outplot
fclose(output_txt);
fileID = fopen(outfile_txt,'r');
formatSpec = '%d %f %d';

A = fscanf(fileID,formatSpec,[3 Inf]);
B=A';
disp(size(B));
figure; hold on;

for bd=2:length(B)
    rangex = [B(bd-1,1) B(bd,1)];
    rangey = [B(bd-1,2) B(bd,2)];
    plot(rangex,rangey,'b.-'), xlim([0 length(B)]), ylim([0 200]);
    if B(bd,3) ~= 0
        plot(B(bd,1),B(bd,2),'ro');
        str = sprintf('(%d, %f)', B(bd,1), B(bd,2));
        text(B(bd,1), B(bd,2)+5, str, 'Color', 'r', 'FontSize', 7);
    end
end
outfile_plot = strcat(dirpath_outplot,outfile);
outfile_plot = strcat(outfile_plot,'.png');
print(outfile_plot,'-dpng')
fclose(fileID);

ss = ss + 1;
if ss > samplesize
    break;
end % if samplesize=5, it will only analyze the first 5 images in each species folder

end % of looping on files in each species folder

end % of looping on all species folders

%%% end of main loop

fclose('all');
