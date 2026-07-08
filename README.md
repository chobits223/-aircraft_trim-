# -aircraft_trim-
这是一个固定翼飞机定直平飞配平计算程序，功能如下：

1.气动建模：基于迎角和侧滑角插值获取升力、阻力、俯仰力矩及舵效导数，支持纵向与侧向气动系数。

2.配平求解：采用牛顿迭代法求解三个平衡方程（升力=重力、俯仰力矩=0、推力=阻力），未知量为飞行速度、升降舵偏角和推力。

3.数值优化：包含雅可比矩阵有限差分、边界限制、线搜索和阻尼因子自适应调整，提高收敛鲁棒性。

4.结果输出：显示配平后的速度、舵偏、推力及迭代历程，并绘制残差范数和状态量的收敛曲线，验证力与力矩平衡。

5.注：所有气动数据与飞机参数均为示例，使用时需要替换为真实值。


This is a trim calculation program for steady straight‑level flight of a fixed‑wing aircraft. Its functions are as follows:

1.Aerodynamic modeling: Interpolates lift, drag, pitching moment, and control effectiveness derivatives based on angle of attack and sideslip angle, supporting both longitudinal and lateral aerodynamic coefficients.

2.Trim solution: Uses the Newton‑Raphson iteration method to solve three equilibrium equations (lift = weight, pitching moment = 0, thrust = drag), with unknowns being flight speed, elevator deflection, and thrust.

3.Numerical optimization: Includes finite‑difference Jacobian computation, boundary constraints, line search, and adaptive damping factor adjustment to improve convergence robustness.

4.Result output: Displays the trimmed speed, control deflection, thrust, and iteration history, and plots convergence curves of residual norm and state variables to verify force and moment balance.

5.Note: All aerodynamic data and aircraft parameters are examples and should be replaced with actual values when used.
