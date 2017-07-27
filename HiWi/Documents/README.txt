

	||-----------------------------------------------------------------------||   
	||-----------------------------------------------------------------------||
	||    Manual for Trommel GUI:					         ||
	||									 ||
	||    Final Update: 24/04/2016						 ||
	||									 ||
	||    Written by: Jiho Yang						 ||
	||                M.Sc. Candidate, Computational Science & Engineering   ||
	||	          Technische Universitat Munchen, Germany	         ||
	||-----------------------------------------------------------------------||
	||-----------------------------------------------------------------------||



사용법 (Damping coefficient vs Velocity):
*실험 데이터와 시뮬레이션 데이터 비교 메뉴얼은 맨 뒤를 참조



1. Browse panel로 불러오고자 하는 데이터 파일 선택 

-File format: .xlsx (다른 파일 포멧도 선택 가능하나 현 버젼에서는 데이터 축출이 불가)

-파일을 불러오면 Z DISPLACEMENT (4번째 센서), Z ACCELERATION (5번, 6번 센서) vs Time (sec) 플롯이 그려짐



2. Offset과 Correction Factor 입력

- Equation used: DEFLECTION = (DISPLACEMENT + OFFSET)*CORRECTION_FACTOR

*주의: 사용된 식을 고려하여 Offset의 +- sign 확인




3. Calibrate panel 클릭 

- DEFLECTION을 계산 후 플롯을 그림

**주의: 다음 단계로 넘어가기 전에 Calibrate panel을 꼭 클릭할 것. Calibration이 되지 않은 상태에서 Calculate Damping panel을 클릭하면 에러 발생. 하지만 혹시 실수로 넘어가서 에러가 발생해도 그냥 Calibrate버튼을 누르면 문제없이 진행됨. Calibrate panel을 여러번 클릭해도 같은 값을 주기때문에 눌렀는지 확신이 들지 않을시 그냥 클릭해도 문제가 안됨



4. Velocity값과 Force (Wheel load)값을 입력

***주의: Velocity와 Force값을 잊지 말고 꼭 입력/변경 할 것. Calculate Damping을 클릭 시 table에 velocity와 force값이 저장되기에 꼭 올바른 값을 입력해 줘야 함. 한번이라도 실수할 시 table을 reset해야할 수도 있음




5. Calculate Damping panel 클릭

- Delta와 damping coefficient 값을 계산 후, delta, velocity, damping coefficient, force값들을 matrix에 저장 후 테이블에 display함




6. 1 - 5번을 주어진 하중값에 대해 반복 (= 하나의 damping coefficient vs velocity 그래프를 그리는데 필요한 데이터가 모일때까지 반복)




7. 주어진 한 하중에 대한 그래프를 그릴 데이타가 다 모이면 Plot Damping panel클릭

- Damping coefficient vs velocity 플롯을 그림

***주의: 주어진 한 하중에 대한 그래프를 그릴 데이타가 다 모이면 꼭 꼭 꼭 plot damping을 클릭해서 플롯을 한번 그릴 것. 그래프를 그리지 않고 다른 하중에 대한 데이타를 부를 시 잘못된 플롯이 그려짐 (=reset을 눌러야 할 수도 있음)

******주의: Damping coefficient vs veloity 플롯이 그려진 figure (figure(6))는 절대 끄지 말 것. 실수로 껐을시 처음부터 다시 반복해야함.



8. 모든 하중에 대한 damping coefficient vs velocity플롯을 생성할 떄까지 반복

*혹 실수로 잘못된 데이터를 테이블에 저장시 Reset버튼 클릭



9. 모든 하중과 속도에 대한 데이터가 모였으면 Export Data으로 데이터 축출(.txt)

***Data column convention: Delta, Velocity, Damping Coefficient
*저장폴더: faml_Trommel_GUI.m이 저장되어 있는 폴더 내에서 Exported Data_MATLAB 폴더 안에 저장됨


%----------------------------------------------------------------------------%


사용법 (Deflection 값에 대한 실험과 시뮬레이션 값 비교):

1. 위 메뉴얼의 1-3번 반복 (calibration까지)



2. 시뮬레이션 값 불러오기



3. Time / Deflection / Scale 중 하나만을 선택 후 Value에 값을 입력하여 Calibrate (Deflection) 클릭

*주의: 셋 중에 하나만을 선택할 것. 혹 실수로 여러개를 선택해도 문제되진 않음 (간략한 메세지만 나타남).
*Net calibration value는 legend에 계속해서 업데이트 되며 표시됨
*Net calibration value 순서 (legend에서): x shift, y shift, y scale








