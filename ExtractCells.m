function cellList = ExtractCells(bwLabel,im,shrinkRate)

im(im>0)=1;
[xdim,ydim]=size(im);
nbr=[1,1;1,0;1,-1;0,1;0,-1;-1,1;-1,0;-1,-1];

epImg = bwmorph(im, 'endpoint');
labelImg=epImg+im;
cellNum = nnz(epImg)/2;
if(mod(nnz(epImg),2)==1)
    disp('error in number of endpoints');
    keyboard
end
clear epImg

ep=find(labelImg==2);
currentCellNum=0;
cellList=cell(1,cellNum);

while(~isempty(ep))
    currentCellNum=currentCellNum+1;
    
    % extract a new end point
    [xt,yt]=ind2sub([xdim ydim],ep(1));
    labelImg(xt,yt)=0;
    
    oep=0; % flag indicating if other endpoint is removed in iteration
    
    % extract the cell body
    tmpPixInd=zeros(50,2);
    pixNum=1;
    tmpPixInd(1,1)=xt; tmpPixInd(1,2)=yt;
    flag=1;
    while(flag)
        flag=0;
        for i=1:8
            x0=xt+nbr(i,1);
            y0=yt+nbr(i,2);
            if(x0<1 || x0>xdim || y0<1 || y0>ydim)
                continue;
            end
            if(labelImg(x0,y0)>0)
                flag=1;
                if(labelImg(x0,y0)==2)
                    if(oep==0)
                        oep=1;
                    else
                        disp('more then two endpoints are removed in one iteration');
                        keyboard
                    end
                end
                
                pixNum=pixNum+1;
                tmpPixInd(pixNum,1)=x0;tmpPixInd(pixNum,2)=y0;
                labelImg(x0,y0)=0;
                xt=x0; yt=y0;
                
                break;
            end
        end
    end
    
    if(oep==0)
        disp('only one endpoint is removed in this iteration');
        keyboard
    end
    
    offset = floor(pixNum*shrinkRate/2)+1;
    cellLength = floor((1-shrinkRate)*pixNum);
    pts=zeros(cellLength,2);
    pts(:,:)=tmpPixInd(offset:1:offset+cellLength-1,:);
    
    thickness = cellThickness(tmpPixInd(1:1:pixNum,:),...
        (bwLabel==bwLabel(pts(1,1),pts(1,2))),xdim,ydim);
    
    cellList{currentCellNum}=struct('pts',pts,'thickness',thickness,'length',0,...
        'targetLength',pixNum,'strip1',[],'strip2',[],'region',[],'intensity',...
        [],'normvec',[],'copyLength',0);
    
    if(~isCloseToBoundary(pts,xdim,ydim))
        cellList{currentCellNum}.copyLength = pixNum;
    end

    ep=find(labelImg==2); 
end

if(currentCellNum~=cellNum)
    disp('error in the number of detected cells');
    disp([currentCellNum,cellNum]);
    keyboard
end
