% u=up line
% l=low line
% b=x axes
% a=alpha for transparency  [0 1]

%varargin: 'significance',pp
       %   'color'
% pp=pval
% CALL IS AS:   colorBetweenLine (freq,coherence,surrogates,alpha,'color',rgb('lightSeaGreen'))


function [] = colorBetweenLine(b,u,l,a,varargin)

if size(u,1)>size(u,2)
   u=u';
end
if size(l,1)>size(l,2)
    l=l';
end
if size(b,1)>size(b,2)
    b=b';
end
 yellow=[1 1 102/255];
 orange=[1 178/255 102/255];
 red= [1 102/255 102/255];
 green=[0 1 0];
 grey=[0.8 0.8 0.8];
 
 if strcmp(varargin{1},'significance')
    pp= varargin{2};
    for i=1:length(pp)-1
        if pp(i)<0.05 &&  pp(i)>=0.01 
           patch(h,[b(i) b(i+1) b(i+1) b(i)], [l(i) l(i+1) u(i+1) u(i)],yellow, 'EdgeColor','none');   
           eval(['alpha ' num2str(a)])
        elseif pp(i)<0.01 &&  pp(i)>=0.001
           patch(h,[b(i) b(i+1) b(i+1) b(i)], [l(i) l(i+1) u(i+1) u(i)],orange , 'EdgeColor','none');   
           eval(['alpha ' num2str(a)])
         elseif pp(i)<0.001
           patch(h,[b(i) b(i+1) b(i+1) b(i)], [l(i) l(i+1) u(i+1) u(i)],red, 'EdgeColor','none');       
           eval(['alpha ' num2str(a)])
        end
    end

 elseif strcmp(varargin{1},'color')
     colin=varargin{2};
       patch(h,[b fliplr(b)], [u fliplr(l)],rgb(colin), 'EdgeColor','none');  
    eval(['alpha ' num2str(a)])
 end

