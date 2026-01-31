/**********************************************************************************************
 *                                                                                            *
 *       CalculiX adapter for heat transfer coupling and mechanical FSI using preCICE         *
 *       Heat transfer adapter developed by Lucía Cheung with the support of SimScale GmbH    *
 *                                                                                            *
 *       Adapter extended to fluid-structure interaction by Alexander Rusch                   *
 *                                                                                            *
 *********************************************************************************************/

#ifndef PRECICEINTERFACE_H
#define PRECICEINTERFACE_H

#include <string.h>
#include "ConfigReader.h"
#include "CCXHelpers.h"

/*
 * PreciceInterface: Structure with all the information of a coupled surface
 * Includes data regarding the surface mesh(es) and the coupling data
 */
typedef struct PreciceInterface {

	char * name;
	int dim;

	// Interface nodes
	int numNodes;
	int * nodeIDs;
	double * nodeCoordinates;
	int nodeSetID;
	int * preciceNodeIDs;
	int nodesMeshID;
	char * nodesMeshName;

	// Interface face elements
	int numElements;
	int * elementIDs;
	int * faceIDs;
	double * faceCenterCoordinates;
	int faceSetID;
	int faceCentersMeshID;
	char * faceCentersMeshName;
	int * preciceFaceCenterIDs;
	int * triangles;

	// Arrays to store the coupling data
	double * nodeScalarData;
	double * nodeVectorData; //Forces, displacements, velocities, positions and displacementDeltas are vector quantities
	double * faceCenterData;

	// preCICE Data IDs
	int temperatureDataID;
	int fluxDataID;
	int kDeltaWriteDataID;
	int kDeltaTemperatureWriteDataID;
	int kDeltaReadDataID;
	int kDeltaTemperatureReadDataID;
	int displacementsDataID; //New data ID for displacements
	int displacementDeltasDataID; //New data ID for displacementDeltas
	int positionsDataID; //New data ID for positions
	int velocitiesDataID; //New data ID for velocities
	int forcesDataID; //New data ID for forces

	// Indices that indicate where to apply the boundary conditions / forces
	int * xloadIndices;
	int * xbounIndices;
	int * xforcIndices;

	// Mapping type if nearest-projection mapping
	int mapNPType;


	int numReadData;
	int numWriteData;
	enum CouplingDataType *readData;
	enum CouplingDataType *writeData;

} PreciceInterface;

/*
 * SimulationData: Structure with all the CalculiX variables
 * that need to be accessed by the adapter in order to do the coupling.
 * A list of variables and their meaning is available in the documentation
 * ccx_2.10.pdf (page 518)
 */
typedef struct SimulationData 
{

	// CalculiX data
	ITG * ialset;						/*---Member of a set or surface. This is a node for node set, an element for element set---*/
	ITG * ielmat;						/*---Contains the material number for element i ---*/
	ITG * istartset;					/*---Pointer to ialset containing the first set number---*/
	ITG * iendset;						/*---Pointer to ialset containing the last set number---*/
	char ** lakon;						/*---Containts label for element i (C3D4, C3D8...)---*/
	ITG * kon;							/*---Containts topology of all elements---*/
	ITG * ipkon;						/*---Points to the location in field kon preceding the topology of element i---*/
	ITG nset;							/*---Number of sets, including surfaces---*/
	char * set;							/*---Name of the set: User defined name---*/
	double * co;						/*---Node coordinates---*/
	ITG nboun;							/*---Total number of boundary conditions (Single Point Constraints)---*/	
	ITG nforc; 							/*---Number of point loads---*/
	ITG * ikboun;						/*---Containts all DOFs of the boundary conditions---*/
	ITG * ikforc; 						/*---Ordered array of the DOFs corresponiding to the point loads---*/
	ITG * ilboun;						/*---Containts all the boundary conditions---*/
	ITG * ilforc; 						/*---Original SPC number for ikforc(i)---*/
	ITG * nelemload;					/*---Element to which distributed load is applied----*/
	int nload;							/*---Number of facial distributed loads---*/
	char * sideload;					/*---Load label; indicated element size to which load is applied---*/
	double nk;							/*---Highest node number---*/
	ITG mt;								/*---Not sure---*/
	double * theta;						/*---Normalized (by tper) size of all previous increments and not including present increment---*/
	double * dtheta;					/*---Normalized (by tper) increment size---*/
	double * tper;						/*---Use given step size---*/
	ITG * nmethod;						/*---Flag that deifnes numerical method: 1: static linear or nonlinear, 2: frequency (linear), 3: buckling, 4: dynamic linear or non-linear, etc---*/
	double * xload;						/*---Concentrated load in direction of idof of node "node" (global coordinates)---*/
	double * xforc; 					/*---Scalar value of the force in one direction---*/
	double * xboun;						/*---Magnitude of constraint at end of a step---*/
	ITG * ntmat_;						/*---Maximum number of temperature data points for any material property for any material---*/
	double * vold;						/*---Displacement of node j in direction i at the start of an iteration---*/
	double * veold;						/*---Velocity of node j in direction i at the start of an iteration---*/
	double * fn;						/*---values of forces read from calculix---*/
	double * cocon;						/*---Conductivity coefficient k at location---*/
	ITG * ncocon;						/*---Number of conductivity constants---*/
	ITG * mi;							/*---Not sure---*/

	// Interfaces
	int numPreciceInterfaces;
	PreciceInterface ** preciceInterfaces;

	// Coupling data
	double * coupling_init_v;
	double coupling_init_theta;
	double coupling_init_dtheta;
	double precice_dt;
	double solver_dt;

} SimulationData;



/**
 * @brief Initialize the preCICE coupling for a CalculiX-based solver.
 *
 * This routine sets up the preCICE adapter and prepares all coupling
 * interfaces required for fluid–structure (or multi-physics) interaction.
 * It performs the following steps:
 *
 *  - Reads the adapter YAML configuration file corresponding to the given
 *    preCICE participant.
 *  - Creates the preCICE solver interface using the preCICE XML configuration.
 *  - Constructs and initializes all coupling interfaces defined in the YAML
 *    file (e.g., surface meshes, read/write data mappings).
 *  - Allocates internal buffers required for coupling state initialization.
 *  - Calls the preCICE initialization routine and retrieves the initial
 *    coupling time step size.
 *  - Initializes coupling data on all registered interfaces.
 *
 * After this function returns, the SimulationData structure is fully
 * configured for use in a preCICE coupling loop (i.e., calls to
 * Precice_ReadCouplingData(), Precice_WriteCouplingData(), and
 * Precice_Advance()).
 *
 * @param[in]  configFilename
 *     Path to the adapter YAML configuration file. This file specifies
 *     the preCICE XML configuration file and the coupling interfaces
 *     associated with the given participant.
 *
 * @param[in]  participantName
 *     Name of the preCICE participant. This must match the participant
 *     name defined in the preCICE XML configuration.
 *
 * @param[in,out] sim
 *     Pointer to the SimulationData structure holding the solver state.
 *     On entry, this structure must contain valid mesh and model data.
 *     On exit, it is populated with all preCICE-related data structures,
 *     including coupling interfaces, buffers, and the initial coupling
 *     time step.
 *
 * @note
 *     This implementation assumes a serial CalculiX participant
 *     (MPI rank = 0, MPI size = 1). Parallel preCICE coupling would
 *     require modifications to the solver interface creation.
 *
 * @warning
 *     All input pointers must be non-NULL. The function will abort
 *     execution if invalid arguments or configuration data are detected.
 */
void Precice_Setup( char * configFilename, char * participantName, SimulationData * sim );

/**
 * @brief Perform the initial preCICE coupling data handshake.
 *
 * This routine initializes the coupling data exchange between the solver
 * and preCICE before entering the main coupling loop. It establishes a
 * consistent initial state for all coupled quantities by performing the
 * following steps:
 *
 *  - Write the solver's initial coupling data (e.g., displacements,
 *    temperatures) to preCICE.
 *  - Signal to preCICE that the initial data exchange phase is complete.
 *  - Read the initial coupling data provided by the coupled participant
 *    (e.g., interface forces or heat fluxes).
 *
 * After this function returns, both the solver and preCICE are synchronized
 * and all coupling fields are initialized consistently, allowing safe
 * entry into the time-stepping or fixed-point coupling loop.
 *
 * @param[in,out] sim
 *     Pointer to the SimulationData structure containing the solver state
 *     and all preCICE coupling information. The structure is updated in
 *     place with the initial coupling data read from preCICE.
 *
 * @note
 *     The order of operations (write → initialize → read) is mandatory
 *     and follows the canonical preCICE initialization protocol.
 *
 * @warning
 *     This function must be called after precicec_initialize() and before
 *     the first call to Precice_Advance(). Skipping this step may lead to
 *     undefined coupling behavior or preCICE runtime errors.
 */
void Precice_InitializeData( SimulationData * sim );


/**
 * @brief Synchronize the solver time step with the preCICE coupling window.
 *
 * This routine adjusts the structural solver time-stepping parameters to
 * ensure consistency with the time step requested by preCICE. The behavior
 * depends on whether the current analysis is steady-state (linear static)
 * or transient (dynamic):
 *
 * - For steady-state (linear static) analyses, the entire solution is
 *   computed in a single coupling step. The normalized solver time is
 *   reset and the solver time step is set equal to the preCICE coupling
 *   time step.
 *
 * - For transient analyses, the solver time increment is constrained to
 *   not exceed the coupling time window provided by preCICE. The normalized
 *   time step used by CalculiX is adjusted accordingly, and the resulting
 *   physical solver time step is synchronized with the preCICE time step.
 *
 * After this function returns, the solver time step stored in
 * sim->solver_dt is consistent with the current preCICE coupling window
 * and can be safely used in the subsequent solver iteration.
 *
 * @param[in,out] sim
 *     Pointer to the SimulationData structure containing the solver and
 *     coupling state. The fields related to time integration (theta, tper,
 *     dtheta, solver_dt) are updated in place.
 *
 * @note
 *     This implementation currently enforces one solver step per preCICE
 *     coupling window. Allowing solver subcycling within a coupling window
 *     would require additional logic.
 *
 * @warning
 *     This function must be called at the beginning of each coupling
 *     iteration, before advancing the solver or exchanging coupling data.
 */
void Precice_AdjustSolverTimestep( SimulationData * sim );


/**
 * @brief Advance the preCICE coupling state by one solver time step.
 *
 * This routine commits the current solver step to preCICE and advances
 * the global coupling time. It informs preCICE of the physical time
 * increment used by the solver, synchronizes with all coupled
 * participants, and retrieves the next coupling time step recommended
 * by preCICE.
 *
 * The returned coupling time step is stored in sim->precice_dt and
 * should be used to adjust the solver time step in the next coupling
 * iteration.
 *
 * @param[in,out] sim
 *     Pointer to the SimulationData structure containing the solver and
 *     coupling state. The field sim->solver_dt must be set prior to this
 *     call. On return, sim->precice_dt is updated with the next coupling
 *     time step provided by preCICE.
 *
 * @note
 *     This function acts as a synchronization barrier: all coupled
 *     solvers must call precicec_advance() before any participant can
 *     proceed to the next coupling iteration.
 *
 * @warning
 *     This function must be called after writing coupling data and
 *     before reading new coupling data for the next iteration. Calling
 *     this function out of sequence may lead to undefined coupling
 *     behavior or preCICE runtime errors.
 */
void Precice_Advance( SimulationData * sim );


/**
 * @brief Query whether the preCICE coupling is still ongoing.
 *
 * This routine checks the internal preCICE coupling state and indicates
 * whether further coupling iterations or time steps are required.
 * It returns true as long as the coupled simulation has not reached its
 * termination condition (e.g., final coupling time or convergence of an
 * implicit coupling scheme).
 *
 * This function is typically used as the condition of the main coupling
 * loop. The solver should continue to exchange data and advance the
 * coupling state while this function returns true.
 *
 * @return
 *     true  if the coupling is still ongoing and further solver iterations
 *           are required.
 *     false if the coupling has finished and the solver may exit the
 *           coupling loop.
 *
 * @note
 *     The termination condition is managed entirely by preCICE. The solver
 *     does not need to track coupling convergence or final time explicitly.
 */
bool Precice_IsCouplingOngoing();


/**
 * @brief Check whether a rollback to a saved iteration checkpoint is required.
 *
 * This routine queries preCICE to determine whether the solver must restore
 * a previously saved iteration checkpoint. Such a rollback is requested by
 * preCICE during implicit coupling schemes when a coupling iteration fails
 * to converge and must be repeated with updated interface data.
 *
 * If this function returns true, the solver is expected to restore its
 * internal state (e.g., displacements, state variables) from the most recent
 * checkpoint before continuing the coupling process.
 *
 * @return
 *     true  if preCICE requests that the solver read and restore an
 *           iteration checkpoint.
 *     false if no rollback is required and the solver may proceed normally.
 *
 * @note
 *     This function only queries the preCICE state. It does not perform the
 *     checkpoint restore itself; the solver must explicitly read the stored
 *     checkpoint data when this function returns true.
 *
 * @warning
 *     This function is relevant only for implicit coupling schemes.
 *     It should be called after Precice_Advance() and before continuing
 *     with the next solver iteration.
 */
bool Precice_IsReadCheckpointRequired();


/**
 * @brief Check whether the solver must write an iteration checkpoint.
 *
 * This routine queries preCICE to determine whether the solver is required
 * to store its current state as an iteration checkpoint. Such checkpoints
 * are used by preCICE during implicit coupling schemes to enable rollback
 * in case a coupling iteration fails to converge.
 *
 * When this function returns true, the solver should save all state
 * information necessary to exactly restore the current iteration, such as
 * displacements, internal state variables, and any other solver-dependent
 * data required for a rollback.
 *
 * @return
 *     true  if preCICE requests that the solver write an iteration
 *           checkpoint.
 *     false if no checkpoint is required at this point.
 *
 * @note
 *     This function only queries the preCICE state and does not perform
 *     the checkpoint write itself. The solver must explicitly write and
 *     later acknowledge the checkpoint when this function returns true.
 *
 * @warning
 *     This function is relevant only for implicit coupling schemes.
 *     It should be called before advancing the coupling state and before
 *     modifying the solver state for the next iteration.
 */
bool Precice_IsWriteCheckpointRequired();


/**
 * @brief Acknowledge completion of an iteration checkpoint restore.
 *
 * This routine notifies preCICE that the solver has successfully restored
 * its state from a previously saved iteration checkpoint. It must be called
 * after preCICE requests a rollback (via the "read-iteration-checkpoint"
 * action) and after the solver has fully restored all required state
 * variables.
 *
 * Calling this function marks the corresponding preCICE action as fulfilled,
 * allowing the coupling algorithm to proceed to the next iteration or time
 * step.
 *
 * @note
 *     This function does not perform the checkpoint restore itself. It only
 *     signals to preCICE that the restore operation has been completed.
 *
 * @warning
 *     This function is required for implicit coupling schemes. Failing to
 *     call it after restoring a checkpoint may cause the coupling procedure
 *     to stall or result in a preCICE runtime error.
 */
void Precice_FulfilledReadCheckpoint();


/**
 * @brief Acknowledge completion of an iteration checkpoint write.
 *
 * This routine notifies preCICE that the solver has successfully written
 * an iteration checkpoint containing its current state. Such checkpoints
 * are requested by preCICE during implicit coupling schemes to enable
 * rollback in case a coupling iteration fails to converge.
 *
 * Calling this function marks the corresponding preCICE action
 * ("write-iteration-checkpoint") as fulfilled, allowing the coupling
 * algorithm to proceed.
 *
 * @note
 *     This function does not perform the checkpoint write itself. The solver
 *     must explicitly save all necessary state information before calling
 *     this routine.
 *
 * @warning
 *     This function is required for implicit coupling schemes. Failing to
 *     call it after writing a checkpoint may cause the coupling procedure
 *     to stall or result in a preCICE runtime error.
 */
void Precice_FulfilledWriteCheckpoint();


/**
 * @brief Restore the solver state from a previously saved iteration checkpoint.
 *
 * This routine performs the actual rollback requested by preCICE during an
 * implicit coupling scheme. It restores the solver state to the values saved
 * at the beginning of the current coupling iteration, allowing the iteration
 * to be repeated after a failed convergence attempt.
 *
 * Specifically, this function:
 *  - Restores the normalized solver time variable (theta)
 *  - Restores the full solution vector containing nodal degrees of freedom
 *
 * The solver time-step size is not restored here, as it is recomputed and
 * synchronized with preCICE in subsequent calls to the time-step adjustment
 * routine.
 *
 * @param[in,out] sim
 *     Pointer to the SimulationData structure containing the stored iteration
 *     checkpoint data and solver state variables.
 *
 * @param[out] v
 *     Pointer to the solver solution vector to be restored (e.g., nodal
 *     displacements). On return, this vector contains the checkpointed
 *     solution state.
 *
 * @note
 *     This function must be called only when preCICE explicitly requests a
 *     rollback (i.e., when the "read-iteration-checkpoint" action is required).
 *
 * @warning
 *     This routine restores only the data explicitly stored in the iteration
 *     checkpoint. For nonlinear analyses, additional internal state variables
 *     may need to be checkpointed and restored to ensure consistency.
 */
void Precice_ReadIterationCheckpoint( SimulationData * sim, double * v );

/**
 * @brief Save the current solver state as an iteration checkpoint.
 *
 * This routine stores a snapshot of the solver state at the beginning of
 * a coupling iteration. The checkpoint enables rollback during implicit
 * preCICE coupling if a coupling iteration fails to converge.
 *
 * The following solver quantities are saved:
 *  - The normalized solver time variable (theta)
 *  - The full solution vector containing nodal degrees of freedom
 *
 * The solver time-step size is intentionally not stored here, as it is
 * recomputed and synchronized with preCICE in subsequent iterations.
 *
 * @param[in,out] sim
 *     Pointer to the SimulationData structure used to store checkpoint data
 *     and solver state variables.
 *
 * @param[in] v
 *     Pointer to the solver solution vector (e.g., nodal displacements) to
 *     be saved in the iteration checkpoint.
 *
 * @note
 *     This function should be called only when preCICE explicitly requests
 *     a checkpoint write (i.e., when the "write-iteration-checkpoint" action
 *     is required).
 *
 * @warning
 *     For nonlinear analyses, additional history-dependent state variables
 *     may need to be included in the checkpoint to ensure a consistent
 *     rollback.
 */
void Precice_WriteIterationCheckpoint( SimulationData * sim, double * v );


/**
 * @brief Read coupling data from preCICE and apply it as solver boundary conditions.
 *
 * This routine retrieves coupling data provided by preCICE and maps it onto
 * the corresponding CalculiX data structures. The data is read for all
 * configured coupling interfaces and for each data type registered as
 * read-data on those interfaces.
 *
 * Depending on the coupling configuration, this function may read and apply:
 *  - Nodal temperatures (Dirichlet thermal boundary conditions)
 *  - Element-face heat fluxes (Neumann thermal boundary conditions)
 *  - Sink temperatures and heat transfer coefficients for convective film
 *    boundary conditions
 *  - Nodal forces (Neumann mechanical boundary conditions)
 *  - Nodal displacements (Dirichlet mechanical boundary conditions)
 *
 * The function first checks whether preCICE has data available for reading.
 * If so, it reads the data in block form using the preCICE C API and applies
 * it to the appropriate CalculiX arrays (e.g., xboun, xforc, xload) via
 * adapter-specific mapping routines.
 *
 * For the FORCES coupling data type, additional post-processing is performed
 * to compute average force components and total force metrics. These values
 * are printed to the console, and a file named "Forces.csv" is written
 * containing per-node force data. This CSV output is intended strictly for
 * diagnostic and post-processing purposes and has no influence on the solver
 * state or the coupling algorithm.
 *
 * @param[in,out] sim
 *     Pointer to the SimulationData structure containing preCICE interface
 *     definitions, coupling metadata, and solver boundary-condition arrays.
 *
 * @note
 *     This function only modifies solver boundary-condition data when
 *     preCICE reports that read-data is available. If no data is available,
 *     the function returns without modifying the solver state.
 *
 * @warning
 *     The diagnostic output written to "Forces.csv" overwrites any existing
 *     file with the same name and should not be used as part of a production
 *     coupling workflow.
 */
void Precice_ReadCouplingData( SimulationData * sim );


/**
 * @brief Write solver coupling data to preCICE.
 *
 * This routine extracts interface-related solution data from the structural
 * solver and sends it to preCICE according to the coupling configuration.
 * Data is written for all registered coupling interfaces and for each data
 * type configured as write-data on those interfaces.
 *
 * Depending on the coupling setup, this function may write:
 *  - Nodal displacements (absolute)
 *  - Nodal displacement increments relative to the current coupling iteration
 *  - Nodal velocities
 *  - Current nodal positions (undeformed coordinates plus displacements)
 *  - Nodal forces
 *
 * The function first checks whether preCICE requires data to be written at
 * the current solver time step or whether an initial data write has been
 * requested. If neither condition is met, no data is written and the
 * function returns immediately.
 *
 * For each data type, solver state variables are gathered from the global
 * CalculiX arrays (e.g., vold, veold, fn, co) and mapped to preCICE interface
 * buffers using adapter-specific helper routines before being sent via the
 * preCICE C API.
 *
 * If preCICE requests the special "write-initial-data" action, this routine
 * fulfills the request after completing the data transfer.
 *
 * @param[in,out] sim
 *     Pointer to the SimulationData structure containing solver state,
 *     preCICE interface definitions, and coupling metadata.
 *
 * @note
 *     The specific quantities written are determined by the adapter
 *     configuration and may differ between coupling interfaces.
 *
 * @warning
 *     This function must be called only after the solver state (displacements,
 *     velocities, forces, etc.) has been updated for the current iteration.
 *     Writing inconsistent or partially updated data may lead to coupling
 *     instability or convergence issues.
 */
void Precice_WriteCouplingData( SimulationData * sim );


/**
 * @brief Release all preCICE coupling resources and finalize the adapter.
 *
 * This routine frees all memory and data structures allocated by the preCICE
 * adapter during setup and coupling. It cleans up iteration checkpoint storage,
 * deallocates all per-interface coupling buffers, and releases the interface
 * objects themselves.
 *
 * After all adapter-owned memory has been freed, the function finalizes the
 * preCICE participant by calling the preCICE C API finalization routine.
 * This signals a clean shutdown of the coupling and releases all internal
 * preCICE communication resources.
 *
 * This function must be called exactly once at the end of the simulation,
 * after all coupling iterations have completed and no further preCICE calls
 * will be made.
 *
 * @param[in,out] sim
 *     Pointer to the SimulationData structure containing preCICE-related
 *     coupling data and interface definitions.
 *
 * @warning
 *     Calling this function while coupling is still ongoing or before all
 *     required preCICE actions have been fulfilled may lead to undefined
 *     behavior or deadlocks.
 *
 * @note
 *     This function frees only adapter-allocated resources. Core solver memory
 *     (e.g., displacement arrays, force arrays) is managed and freed separately
 *     by the solver.
 */
void Precice_FreeData( SimulationData * sim );


/**
 * @brief Initialize and configure a preCICE coupling interface.
 *
 * This routine sets up a single coupling interface between CalculiX and
 * preCICE based on the interface configuration provided via the adapter
 * configuration file. It initializes all interface data structures, configures
 * the required preCICE meshes (nodal and/or face-center meshes), and registers
 * all coupling data to be read from or written to preCICE.
 *
 * The setup process includes:
 *  - Querying the spatial dimension of the coupling from preCICE
 *  - Initializing all interface pointers and preCICE data identifiers
 *  - Storing the interface patch name and mapping type
 *  - Configuring the nodal mesh used for displacement-, force-, or
 *    velocity-based coupling
 *  - Optionally configuring a face-center mesh and associated surface
 *    triangulation for face-based coupling quantities (e.g., thermal fluxes)
 *  - Registering read and write coupling data and obtaining preCICE data IDs
 *
 * After successful completion, the interface is fully configured and ready
 * for data exchange during the coupling loop.
 *
 * @param[out] interface
 *     Pointer to an allocated PreciceInterface structure that will be
 *     initialized and populated by this function.
 *
 * @param[in] sim
 *     Pointer to the SimulationData structure containing solver state and
 *     CalculiX model information required to build the coupling meshes.
 *
 * @param[in] config
 *     Pointer to the InterfaceConfig structure defining the coupling interface,
 *     including patch name, mesh names, mapping type, and read/write data
 *     specifications.
 *
 * @note
 *     The nodal and face-center meshes are configured only if corresponding
 *     mesh names are provided in the interface configuration.
 *
 * @warning
 *     The PreciceInterface structure must be freed using
 *     PreciceInterface_FreeData() after coupling has finished to avoid memory
 *     leaks.
 */
void PreciceInterface_Create(
    PreciceInterface * interface,
    SimulationData * sim,
    InterfaceConfig const * config
);


/**
 * @brief Configure and register the face-center mesh for a preCICE interface.
 *
 * This function constructs a face-center-based coupling mesh for the given
 * interface using a surface set defined in the CalculiX model. The surface
 * set is derived from the interface patch name and is expected to contain
 * element faces representing the coupling boundary.
 *
 * The routine performs the following steps:
 *  - Determines the CalculiX surface set associated with the interface
 *  - Extracts the element IDs and corresponding face IDs for all faces in
 *    the surface set
 *  - Computes the geometric center of each face using element connectivity
 *    and nodal coordinates
 *  - Registers the face-center coordinates as vertices of a preCICE mesh
 *  - Stores the resulting preCICE vertex IDs for later data exchange
 *
 * The face-center mesh is typically used for face-based coupling quantities
 * such as heat fluxes, sink temperatures, or heat transfer coefficients.
 *
 * @param[in,out] interface
 *     Pointer to the PreciceInterface structure being configured. On return,
 *     this structure contains the face-center mesh ID, face-center coordinates,
 *     element/face mappings, and preCICE vertex IDs.
 *
 * @param[in] sim
 *     Pointer to the SimulationData structure containing CalculiX mesh data,
 *     including connectivity, coordinates, and set definitions.
 *
 * @note
 *     This function assumes tetrahedral elements when computing face centers.
 *     The face-center mesh is configured only if a face-center mesh name was
 *     provided in the interface configuration.
 *
 * @warning
 *     All memory allocated within this function must be released via
 *     PreciceInterface_FreeData() to avoid memory leaks.
 */
void PreciceInterface_ConfigureFaceCentersMesh(
    PreciceInterface * interface,
    SimulationData * sim
);


/**
 * @brief Configure and register the nodal coupling mesh for a preCICE interface.
 *
 * This function constructs the nodal coupling mesh associated with a given
 * interface by identifying the corresponding CalculiX node set, extracting
 * the nodal coordinates, and registering these nodes as vertices in preCICE.
 * The resulting mesh is used for node-based coupling quantities such as
 * displacements, displacement deltas, velocities, and nodal forces.
 *
 * The setup process includes:
 *  - Determining the CalculiX node set associated with the interface patch name
 *  - Extracting the node IDs belonging to this set
 *  - Computing the current nodal coordinates (including deformation, if any)
 *  - Registering the nodal coordinates as vertices of a preCICE mesh
 *  - Storing the resulting preCICE vertex IDs for later data exchange
 *  - Optionally configuring nodal connectivity when nearest-projection mapping
 *    is requested
 *
 * @param[in,out] interface
 *     Pointer to the PreciceInterface structure being configured. On return,
 *     this structure contains the nodal mesh ID, nodal coordinates, node IDs,
 *     and corresponding preCICE vertex IDs.
 *
 * @param[in] sim
 *     Pointer to the SimulationData structure containing CalculiX mesh data,
 *     including nodal coordinates, displacements, set definitions, and
 *     connectivity information.
 *
 * @note
 *     The node IDs are referenced directly from the global CalculiX set array
 *     (ialset). No deep copy of the node list is performed.
 *
 * @note
 *     Nearest-projection mapping requires additional connectivity information
 *     and is enabled only if mapNPType is set to 1 in the interface
 *     configuration.
 *
 * @warning
 *     Memory allocated within this function must be released via
 *     PreciceInterface_FreeData() to avoid memory leaks.
 */
void PreciceInterface_ConfigureNodesMesh(
    PreciceInterface * interface,
    SimulationData * sim
);


/**
 * @brief Terminate execution if the nodes mesh ID is not valid
 * @param interface
 */
void PreciceInterface_EnsureValidNodesMeshID( PreciceInterface * interface );

/**
 * @brief Configures the faces mesh (for tetrahedral elements only)
 * @param interface
 * @param sim
 */
void PreciceInterface_ConfigureTetraFaces( PreciceInterface * interface, SimulationData * sim );

/**
 * @brief Construct surface connectivity for nearest-projection mapping in preCICE.
 *
 * This function builds the geometric surface connectivity required by preCICE
 * when nearest-projection (NP) mapping is enabled for an interface. It derives
 * the interface surface from a CalculiX face (surface) set, extracts all
 * surface faces, and reconstructs the triangular surface mesh by identifying
 * the nodal connectivity of each face.
 *
 * The routine performs the following steps:
 *  - Identifies the CalculiX surface set associated with the interface patch
 *  - Determines the number of surface faces in the set
 *  - Extracts the (element ID, face ID) pairs defining the surface
 *  - Reconstructs triangular face connectivity for each surface face
 *  - Stores the resulting triangle node indices for use by preCICE
 *
 * The resulting triangulation is used by preCICE to perform geometric
 * nearest-projection mapping of coupling quantities (e.g., forces,
 * displacements) between solver interfaces.
 *
 * @param[in,out] interface
 *     Pointer to the PreciceInterface structure being configured. On return,
 *     this structure contains the surface triangulation, face IDs, and element
 *     IDs required for nearest-projection mapping.
 *
 * @param[in] sim
 *     Pointer to the SimulationData structure containing CalculiX mesh data,
 *     including element connectivity, node coordinates, and surface set
 *     definitions.
 *
 * @note
 *     This function assumes tetrahedral elements, for which each surface face
 *     is a triangle. Support for other element types would require extending
 *     the face reconstruction logic.
 *
 * @note
 *     This routine is called only when nearest-projection mapping is enabled
 *     (mapNPType == 1).
 *
 * @warning
 *     Memory allocated within this function must be released via
 *     PreciceInterface_FreeData() to avoid memory leaks.
 */
void PreciceInterface_NodeConnectivity(
    PreciceInterface * interface,
    SimulationData * sim
);


/**
 * @brief Configure read/write coupling data for a preCICE interface.
 *
 * This function initializes all coupling data channels associated with a given
 * interface based on the adapter configuration. It allocates data buffers,
 * resolves preCICE data IDs, and constructs index mappings that link preCICE
 * coupling data to the appropriate CalculiX data structures.
 *
 * The configuration supports both node-based and face-based coupling data,
 * including (but not limited to):
 *  - Nodal data: Temperature, Forces, Displacements, Displacement Deltas,
 *    Velocities, Positions
 *  - Face-based data: Heat Flux, Sink Temperature (film BC), Heat Transfer
 *    Coefficient
 *
 * For each requested read or write quantity, this routine:
 *  - Determines whether the data is node-based or face-based
 *  - Ensures the required preCICE mesh (nodes or face centers) is available
 *  - Retrieves the corresponding preCICE DataID
 *  - Allocates temporary buffers for data exchange
 *  - Builds index maps (xbounIndices, xforcIndices, xloadIndices) that connect
 *    preCICE data entries to CalculiX arrays such as xboun, xforc, and xload
 *
 * After successful execution, the PreciceInterface object contains all
 * information required to perform data exchange during the coupling loop via
 * Precice_ReadCouplingData() and Precice_WriteCouplingData().
 *
 * @param[in,out] interface
 *     Pointer to the PreciceInterface structure being configured. On return,
 *     this structure contains:
 *       - Registered preCICE DataIDs for all coupling quantities
 *       - Allocated buffers for scalar and vector data
 *       - Index mappings linking preCICE data to CalculiX arrays
 *       - Enumerations identifying read/write data types
 *
 * @param[in] sim
 *     Pointer to the SimulationData structure containing CalculiX model data,
 *     including loads, boundary conditions, node and element sets, and solver
 *     arrays required for index mapping.
 *
 * @param[in] config
 *     Pointer to the InterfaceConfig structure defining the coupling data for
 *     this interface, including the names and number of read and write data
 *     channels as specified in the adapter configuration file.
 *
 * @note
 *     Node-based data requires a valid nodal preCICE mesh, while face-based data
 *     requires a face-center mesh. These meshes must be configured prior to
 *     calling this function.
 *
 * @warning
 *     All memory allocated within this routine must be released via
 *     PreciceInterface_FreeData() to avoid memory leaks.
 *
 * @warning
 *     Any unknown or unsupported read/write data name in the configuration will
 *     cause the program to terminate with an error.
 */
void PreciceInterface_ConfigureCouplingData(
    PreciceInterface * interface,
    SimulationData * sim,
    InterfaceConfig const * config
);


/**
 * @brief Frees the memory
 * @param preciceInterface
 */
void PreciceInterface_FreeData( PreciceInterface * preciceInterface );

void PreciceInterface_CollectData( PreciceInterface * preciceInterface );


#endif // PRECICEINTERFACE_H
