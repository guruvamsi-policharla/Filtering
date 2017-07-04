function clear_pane_axes(h)
%Deletes all axes present on a pane
child_handles = allchild(h);
    for i = 1:size(child_handles,1)
        if(strcmp(get(child_handles(i),'Type'),'axes'))
            delete(child_handles(i))
        end
    end