function [ImOut, nrObjects] = iCantSeeSh_t(Im)
% This code was written by Simon Aertssen s181603 for the course 02502: 
% Image Analysis. The user 

    %% Setup
    % Rescale for computational purposes.
    [R, K] = size(Im);

    Im = imresize(Im,[NaN, 1000]); 
    [r, k, l] = size(Im);
    
    %% Colour analysis: binary image
    % According to the DTU style guide, sign colour is HEX = 990000, HSV =
    % (0,100,60).
    
    ImHSV = rgb2hsv(Im);
    H = ImHSV(:, :, 1); S = ImHSV(:, :, 2); V = ImHSV(:, :, 3);

    ImThresh = (H < 0.03 | H > 0.85) & S > 0.45;
    ImWhite = ( (H < 0.2 | H > 0.9) & S < 0.5 & V > 0.55 ) ;
    
    %% Morphology: closing the image
    % Get rid of noise: open the image with a small st element. We wish to
    % keep elements that are close together separate, so no dilation yet.
    
    st = 2;
    box = strel('square', st);
    ImRed = imopen(ImRed, box);

    %% Label groups
    [ImLab, LabCount] = bwlabel(ImRed, 8);
    
    %% Is there a possible sign? Filter on properties
    % Gather properties in readable format
    Props = regionprops(ImLab, 'Area', 'BoundingBox', 'Image');
    
    
    % Prepare ImLab to be used as control map
    box = strel('square', 4*st); Shapes = imclose(ImLab, box);
    edge = strel('square', 3*st); Shapes = imerode(Shapes, edge);
    

    for i = 1:LabCount
        
        % Gather all pixels for form i
        n = find(ImLab == i);

        % Sort on areas:
        Area = Props(i).Area;
        if Area/(r*k) > 3.3333e-04 

           % If not too small, see if bounding box has right dimensions:
           BBox = round( vec2mat(Props(i).BoundingBox, 4) );
           if (BBox(3)/BBox(4) > 2) & (BBox(3)/BBox(4) < 6.9)

              % Now see if there are white letters inside of the form of
              % our possible sign: two conditions!
                whitePxls = 0;
                
                for x = BBox(2) + (0:BBox(4)) 
                    for y = BBox(1) + (0:BBox(3))
                        if (Shapes(x, y) == i) && (ImWhite(x , y ) == 1) 
                             whitePxls = whitePxls + 1;
                        end
                    end
                end
                
                if whitePxls/Area > 0.01
                    
                else 
                    LabCount = LabCount - 1;
                    ImLab(n) = 0; 
                    continue
                end
                
           else 
                LabCount = LabCount - 1;
                ImLab(n) = 0;
                continue
           end

        else
            LabCount = LabCount - 1;
            ImLab(n) = 0; 

        end
    end
    
    if LabCount == 1
        fprintf('One sign has been detected \n')
        
    else
        fprintf('%d signs have been detected \n', LabCount)
    end
    
     %% Rescale to obtain original format
     ImLab = imresize(ImLab,[R, NaN]); 
    

    %% Apply ConvexHull for smooth edges
    CH = bwconvhull(ImLab,'objects');
    ImLab = imclose(CH, box);
    
    %Relabel to obtain 0, 1, 2, ...
    ImLab = bwlabel(ImLab, 8);
    
    %% Output
    ImOut = ImLab;
    nrObjects = LabCount;
    
end
