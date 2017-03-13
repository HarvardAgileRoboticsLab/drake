if(lcm_FOUND)
  # Defines a shared library for use by KUKA iiwa demos based on System 2.0.
  add_library_with_exports(LIB_NAME drakeKukaIiwaArmCommon SOURCE_FILES
      iiwa_common.cc iiwa_lcm.cc)
  target_link_libraries(drakeKukaIiwaArmCommon
    drakeIK
    drakeLCMSystem2
    drakeRBM
    drakeMultibodyParsers)
  drake_install_libraries(drakeKukaIiwaArmCommon)
  drake_install_headers(
      iiwa_common.h iiwa_lcm.h)
  drake_install_pkg_config_file(drake-kuka-iiwa-arm-common
      TARGET drakeKukaIiwaArmCommon
      LIBS -ldrakeKukaIiwaArmCommon
      REQUIRES drake-lcm-system drake-lcmtypes-cpp drake-rbm)

  drake_install_headers(kuka_message_handler.h)

  drake_install_headers(kuka_lcm_planner.h)
  drake_install_headers(parse_json.h)
  drake_install_headers(rigid_body_plan_parser.h)
  #add_executable(test_parse_json test_parse_json.cc)
  add_executable(run_kuka_dircol_planner run_kuka_dircol_planner.cc)
  target_link_libraries(run_kuka_dircol_planner
      drakeKukaIiwaArmCommon
      drakeTrajectories
      drakeMultibodyParsers
      drakeIK
      drakeRBM
      drakeLCMTypes
      lcm
      robotlocomotion-lcmtypes-cpp)

  add_executable(run_kuka_iiwa_arm_dynamics run_kuka_iiwa_arm_dynamics.cc)
  target_link_libraries(run_kuka_iiwa_arm_dynamics
      drakeKukaIiwaArmCommon
      drakeMultibodyParsers
      drakeRigidBodyPlant
      drakeSystemAnalysis
      drakeSystemControllers
      drakeLCMSystem2
      gflags)
  drake_add_test(
      NAME run_kuka_iiwa_arm_dynamics
      COMMAND run_kuka_iiwa_arm_dynamics --duration 0.5
      SIZE medium)

  add_executable(kuka_ik_demo kuka_ik_demo.cc)
  target_link_libraries(kuka_ik_demo
    drakeIK
    drakeMultibodyParsers
    drakeRBM
    robotlocomotion-lcmtypes-cpp
    lcm)

  add_executable(kuka_plan_runner kuka_plan_runner.cc)
  target_link_libraries(kuka_plan_runner
    drakeKukaIiwaArmCommon
    drakeTrajectories
    drakeMultibodyParsers
    drakeRBM
    drakeLCMTypes
    robotlocomotion-lcmtypes-cpp
    lcm)
    
  add_executable(kuka_plan_runner_draper kuka_plan_runner_draper.cc)
  target_link_libraries(kuka_plan_runner_draper
    drakeKukaIiwaArmCommon
    drakeTrajectories
    drakeRBM
    drakeLCMTypes
    robotlocomotion-lcmtypes-cpp
    lcm)



  add_executable(kuka_plan_runner_inertia_J5toJ7_chirp kuka_plan_runner_inertia_J5toJ7_chirp.cc)
  target_link_libraries(kuka_plan_runner_inertia_J5toJ7_chirp
    drakeKukaIiwaArmCommon
    drakeTrajectories
    drakeMultibodyParsers
    drakeRBM
    drakeLCMTypes
    robotlocomotion-lcmtypes-cpp
    lcm)

  add_executable(kuka_plan_runner_inertia_J2toJ7_chirp kuka_plan_runner_inertia_J2toJ7_chirp.cc)
  target_link_libraries(kuka_plan_runner_inertia_J2toJ7_chirp
    drakeKukaIiwaArmCommon
    drakeTrajectories
    drakeMultibodyParsers
    drakeRBM
    drakeLCMTypes
    robotlocomotion-lcmtypes-cpp
    lcm)

  add_executable(kuka_plan_runner_friction_estimation kuka_plan_runner_friction_estimation.cc)
  target_link_libraries(kuka_plan_runner_friction_estimation
    drakeKukaIiwaArmCommon
    drakeTrajectories
    drakeMultibodyParsers
    drakeRBM
    drakeLCMTypes
    robotlocomotion-lcmtypes-cpp
    lcm)

  add_executable(kuka_plan_runner_mass_and_CoM_estimation_J5toJ7 kuka_plan_runner_mass_and_CoM_estimation_J5toJ7.cc)
  target_link_libraries(kuka_plan_runner_mass_and_CoM_estimation_J5toJ7
    drakeKukaIiwaArmCommon
    drakeTrajectories
    drakeMultibodyParsers
    drakeRBM
    drakeLCMTypes
    robotlocomotion-lcmtypes-cpp
    lcm)

  add_executable(kuka_plan_runner_mass_and_CoM_estimation_J2toJ7 kuka_plan_runner_mass_and_CoM_estimation_J2toJ7.cc)
  target_link_libraries(kuka_plan_runner_mass_and_CoM_estimation_J2toJ7
    drakeKukaIiwaArmCommon
    drakeTrajectories
    drakeMultibodyParsers
    drakeRBM
    drakeLCMTypes
    robotlocomotion-lcmtypes-cpp
    lcm)
            
  add_executable(kuka_control kuka_control.cc)
  target_link_libraries(kuka_control
    drakeKukaIiwaArmCommon
    drakeTrajectories
    drakeRBM
    drakeLCMTypes
    robotlocomotion-lcmtypes-cpp
    lcm)

  add_executable(kuka_control_gravity_comp_test kuka_control_gravity_comp_test.cc)
  target_link_libraries(kuka_control_gravity_comp_test
    drakeKukaIiwaArmCommon
    drakeTrajectories
    drakeRBM
    drakeLCMTypes
    robotlocomotion-lcmtypes-cpp
    lcm)

  add_executable(kuka_control_torque_comp_test kuka_control_torque_comp_test.cc)
  target_link_libraries(kuka_control_torque_comp_test
    drakeKukaIiwaArmCommon
    drakeTrajectories
    drakeRBM
    drakeLCMTypes
    robotlocomotion-lcmtypes-cpp
    lcm)

  add_executable(kuka_simulation kuka_simulation.cc)
  target_link_libraries(kuka_simulation
    drakeKukaIiwaArmCommon
    drakeMultibodyParsers
    drakeRigidBodyPlant
    drakeSystemAnalysis
    drakeSystemControllers
    drakeLCMSystem2)
  drake_add_test(
      NAME kuka_simulation COMMAND kuka_simulation -simulation_sec 0.1)

#  add_executable(kuka_ik_sine_demo kuka_ik_sine_demo.cc)
#  target_link_libraries(kuka_ik_sine_demo drakeIK
#      drakeTrajectories
#      drakeRBM
#      drakeLCMSystem)

#  add_executable(kuka_id_gravity_comp kuka_id_gravity_comp.cc)
#  target_link_libraries(kuka_id_gravity_comp
#      drakeTrajectories
#      drakeRBM
#      drakeLCMSystem)

#  add_executable(run_kuka_iiwa_gravity_compensated_torque_control
#      run_kuka_iiwa_gravity_compensated_torque_control.cc)
#  target_link_libraries(run_kuka_iiwa_gravity_compensated_torque_control
#      drakeKukaIiwaArm
#      drakeRBSystem
#      drakeLCMSystem
#      gflags)

#  add_executable(run_kuka_iiwa_gravity_compensated_position_control
#      run_kuka_iiwa_gravity_compensated_position_control.cc)
#  target_link_libraries(run_kuka_iiwa_gravity_compensated_position_control
#      drakeKukaIiwaArm
#      drakeRBSystem
#      drakeLCMSystem
#      gflags)

  add_executable(iiwa_wsg_simulation iiwa_wsg_simulation.cc)
  target_link_libraries(iiwa_wsg_simulation
    drakeKukaIiwaArmCommon
    drakeKukaIiwaWorld
    drakeLCMSystem2
    drakeRigidBodyPlant
    drakeSchunkWsg
    drakeSystemAnalysis
    drakeSystemControllers)
  drake_add_test(
      NAME iiwa_wsg_simulation COMMAND iiwa_wsg_simulation -simulation_sec 0.01)

endif()

add_subdirectory(controlled_kuka)
add_subdirectory(iiwa_world)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
