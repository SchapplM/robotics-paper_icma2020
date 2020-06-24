% Return end effector position of the robot
% Input:
% l [6x1]
%   Kinematic parameters
% q [Nx6]
%   Motor positions for N different robot poses
% 
% Output:
% freturn [Nx3]
%   End effector positions for the N input poses of q

% Bjoern Volkmann, bjoern.volkmann@imes.uni-hannover.de, 2020-06
% (C) Institut fuer Mechatronische Systeme, Leibniz Universitaet Hannover

function freturn = posEE(l, q)

	t3866 = sin(q(:,2));
	t3869 = cos(q(:,2));
	t3874 = sin(q(:,5)) .* l(6);
	t3861 = cos(q(:,4)) .* t3874 - l(4);
	t3862 = cos(q(:,5)) .* l(6) + l(5);
	t3865 = sin(q(:,3));
	t3868 = cos(q(:,3));
	t3871 = t3861 .* t3868 + t3862 .* t3865 - l(3);
	t3872 = t3865 .* t3861 - t3862 .* t3868;
	t3876 = t3871 .* t3866 + t3872 .* t3869 - l(2);
	t3873 = sin(q(:,4)) .* t3874;
	t3870 = cos(q(:,1));
	t3867 = sin(q(:,1));
	freturn = [t3867 .* t3873 - t3876 .* t3870 t3876 .* t3867 + t3870 .* t3873 t3872 .* t3866 - t3869 .* t3871 + l(1)];

end
