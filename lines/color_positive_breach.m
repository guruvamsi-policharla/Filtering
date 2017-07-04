function [] = color_positive_breach(h,line_base,line_up, line_down, varargin)
%COLOR_POSITIVE_OVERFLOW Colors the breach above a particular line
%   Plots polygons to color the positive breach in a plot
%   use flipped if you're flipping the axes
%   colors as found in RGB
color = 'black'; 

if(strcmp(varargin{1},'color'))
     color=varargin{2};
end

if(strcmp(varargin{3},'flipped'))
    flipped = 1;
else 
    flipped = 0;
end

line_up = line_up(:);
line_down = line_down(:);
line_base = line_base(:);

L = min(length(line_up),length(line_down));

i = 1;
j = 1;
while i < L
    if line_up(i)>line_down(i)
        while(line_up(i)>=line_down(i) && i<L)
            x(j) = line_base(i);
            y_up(j) = line_up(i);
            y_down(j) = line_down(i);            
            i = i+1;
            j = j+1;
        end        
        y = [y_up,fliplr(y_down)];
        x = [x,fliplr(x)];
        
        if flipped == 0
            patch(h,x,y,rgb(color),'FaceAlpha',.5);
        else
            patch(h,y,x,rgb(color),'FaceAlpha',.5);
        end
        x = [];
        y_up = [];
        y_down = [];
        j = 1;
        
    end
    i = i + 1;
end

