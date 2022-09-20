width=432; %comes from thresholded kymograph in FIJI
height=330; %comes from thresholded kymograph in FIJI

steps=3;
%border=[xborder,yborder,slice];
k=1;
sumset=zeros(10000,6);
cutoff=15000;
for sn=0:19
    tic;
    bordern=border(border(:,3)==sn,1:2);
     xi0=bordern(1,1);
     start=floor((max(bordern,[],1)-min(bordern,[],1))/steps);
     xist=start(1);
     yist=start(2);
     for yi0s = 1:floor(yist/3)
            for xis = floor(xist/4):floor(xist/2)
                for yis = floor(yist/4):floor(yist/2)
                    yi0=steps*yi0s;
                    xi=steps*xis;
                    yi=steps*yis;
                    xline = [ones(1,yi0)*xi0,xi0+1:xi0+xi,ones(1,yi)*(xi+xi0),xi+xi0+1:xi+xi0+xi,ones(1,yi)*(2*xi+xi0),2*xi+xi0+1:2*xi+xi0+max(xi,width-2*xi-xi0)]';
                    yline = [1:yi0,ones(1,xi)*yi0,yi0+1:yi0+yi,ones(1,xi)*(yi+yi0),yi+yi0+1:yi+yi0+yi,ones(1,max(xi,width-2*xi-xi0))*(2*yi+yi0)]';
                    
                    n=1;
                    while ((n<length(xline)-1)*(xline(n)<width)*(yline(n)<height))
                        n=n+1;
                    end
                    
                    xline=xline(1:n);
                    yline=yline(1:n);
                    sumsq=0;
                    for r=1:n-1
                        xpoint=bordern(bordern(:,2)==yline(r),1);
                        ypoint=bordern(bordern(:,1)==xline(r),2);
                        avex=xline(r);
                        avey=yline(r);
                        if ~isempty(xpoint)
                            avex=mean(xpoint);
                        end
                        if ~isempty(ypoint)
                            avey=mean(ypoint);
                        end
                        if xline(r+1)-xline(r)==1
                            sumsq=sumsq+(yline(r)-avey)^2;
                        else
                            sumsq=sumsq+(xline(r)-avex)^2;
                        end
                    end
                    if sumsq<cutoff*10
                        sumset(k,:)=[sumsq/cutoff sn xi0 yi0 xi yi];
                        k=k+1;
                    end
                end
            end
     end
     
            toc
end
sumset=sumset(sumset(:,1)~=0,:);
