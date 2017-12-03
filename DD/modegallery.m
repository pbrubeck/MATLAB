function [] = modegallery(xx,yy,uuu)
nmodes=size(uuu,3);
mode=1;
handle=surf(xx, yy, uuu(:,:,mode));
fprev=get(gcf, 'WindowKeyPressFcn');
function []=KeyboardPressed(obj, evnt)
    switch evnt.Key
        case 'rightarrow'
            mode=mod(mode, nmodes)+1;
        case 'leftarrow'
            mode=mod(mode-2, nmodes)+1;
    end
    set(handle,'ZData',uuu(:,:,mode));
    try
      fprev(obj, evnt);
    catch
    end
end
set(gcf, 'WindowKeyPressFcn', @KeyboardPressed);
end