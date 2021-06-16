function [varargout] = indexEqns(crysName)

%%% Keeping some text around for biaxial crystals
% if nargout == 2
    syms nE(l) nO(l) nE_Theta(l,theta) w
% elseif nargout == 3
%     syms nX(l) nY(l) nZ(l) nE_Theta_Phi(l,theta,phi) w
% else
%     error('Wrong number of outputs');
% end

switch crysName
    case 'BBO'
        nO(l) = sqrt( 2.7359 + 0.01878/((l/10^-6)^2 - 0.01822)...
            - 0.01354 * (l/10^-6)^2 );
        nE(l) = sqrt( 2.3753 + 0.01224/((l/10^-6)^2 - 0.01667)...
            - 0.01516 * (l/10^-6)^2 );
        dNL = 2.01 * 10^-12;
    otherwise
        error('Crystal name not found.');
end
    
nE_Theta(l,theta) = sqrt( 1/ (cosd(theta)^2/nO(l)^2 + sind(theta)^2/nE(l)^2) );

%%% Keeping some text around for biaxial crystals
% if nargout == 3
    varargout{1} = nO;
    varargout{2} = nE_Theta;
    varargout{3} = dNL;
% elseif nargout == 4
%     varargout{1} = nx;
%     varargout{2} = ny;
%     varargout{3} = nz;
%     varargout{4} = dNL;
% else
%     error('Wrong number of outputs');
% end


end