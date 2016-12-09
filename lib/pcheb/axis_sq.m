function axis_sq
%AXIS_SQ  Square axes with the same tick-labels. 

    axis square

    h = gca;  tick = get(h,'XTick');  set(h,'YTick',tick);
    
    xlabel x,  ylabel y

end
