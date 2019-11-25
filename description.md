# Linear Program Solver 
### by Y. Cui, K. Morikuni, T. Tsuchiya and K. Hayami.
First version: Aug 2015 \
Latest update: Nov 2019 

**_This project is licensed under the terms of the GNU license <http://www.gnu.org/copyleft/gpl.html>._**

## Reference
Please find the full article: \
_Implementation of interior-point methods for LP based on 
Krylov subspace iterative solvers with inner-iteration preconditioning_ \
Cui, Y., Morikuni, K., Tsuchiya, T., Hayami, K. \
Comput Optim Appl (2019) 74: 143. \
<https://doi.org/10.1007/s10589-019-00103-y>\
If you use these codes in research for publication, please cite the paper.

## Contents
The repository contains:
* _interior_point_solver_
	* `main.m`: the entry point of solving the interior-point problem
	* `NIILP.m`: the functions of the linear program solver
* _linear_solvers_
	* `ABNESOR4IP_scale.c`: the AB-GMRES solver preconditioned by NE-SOR inner iterations with row scaling.
	* `CGNE4IP_scale.c`: the CGNE solver preconditioned by NE-SSOR inner iterations with row scaling.
	* `MRNE4IP_scale.c`: the MRNE solver preconditioned by NE-SSOR inner iterations with row scaling.\
	To use these C-code inside the interior-point solver `NIILP.m`, please compile them as MEX (MATLAB Executable) files. For more instructions, see <https://uk.mathworks.com/help/matlab/ref/mex.html>.

	The full article and code for AB-GMRES solver with NE-SOR and BA-GMRES solver with NR-SOR can be found in Dr. Keiichi Morikuni's homepage: <https://researchmap.jp/KeiichiMorikuni/Implementations/>. Please cite the relevant paper if you use these linear solvers for publication.
* _html_
	* `main.html`: the html file published by `main.m`
	* `NIILP.html`: the html file published by `NIILP.m`\
	The code in the html files can be extracted using MATLAB function `grabcode url`.

## Contact
Please feel free to get in touch if you have any questions or suggestions: cuiyiran2012@hotmail.com 