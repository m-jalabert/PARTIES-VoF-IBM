function [] = hdf5_reset_to_zero(filename)
% HDF5_RESET_TO_ZERO   Sets time and iteration count to zero for HDF5 file
%    HDF5_RESET_TO_ZERO(filename) Edits filename, setting '/time' to zero
%    and, if it is a resume file, '/noutput' and '/ntime' to zero as well

% For all files
h5write(filename, '/time', 0)

% For Resume_*.h5 files
if (~isempty(strfind(filename,'Resume')))
	h5write(filename, '/noutput', int32(0))
	h5write(filename, '/ntime', int32(0))
	t = h5read(filename, '/timer');
	h5write(filename, '/timer', zeros(size(t)))
end