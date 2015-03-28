/*------------------------------------------------------------------------
 *   initialisation of repeated comunications. This may reduce the
 *   network overhead. Communication is started later
 *   with MPI_START(request)
 *   last update 31/07/00, T. Bohlen
 *  ----------------------------------------------------------------------*/

#include "fd.h"

void comm_ini(float ** bufferlef_to_rig, float ** bufferrig_to_lef, 
float ** buffertop_to_bot, float ** bufferbot_to_top, 
MPI_Request *req_send, MPI_Request *req_rec){


	extern int NX, NY, INDEX[5], FDORDER;
	extern const int TAG1,TAG2,TAG3,TAG4;
	int fdo2;

	

	/* comunication initialisation for persistent communication */

	
	fdo2 = 2*(FDORDER/2 + 1);
	

	/* buffer arrays are copied into local buffers using buffered send (bsend),
	  Actually send is activated by MPI_START(request) within time loop.
	  MPI_BSEND and MPI_RECV (see below) are then non-blocking.
	*/
	MPI_Bsend_init(&bufferlef_to_rig[1][1],NY*fdo2,MPI_FLOAT,INDEX[1],TAG1,MPI_COMM_WORLD,&req_send[0]);
	MPI_Bsend_init(&bufferrig_to_lef[1][1],NY*fdo2-1,MPI_FLOAT,INDEX[2],TAG2,MPI_COMM_WORLD,&req_send[1]);
	MPI_Bsend_init(&buffertop_to_bot[1][1],NX*fdo2,MPI_FLOAT,INDEX[3],TAG3,MPI_COMM_WORLD,&req_send[2]);
	MPI_Bsend_init(&bufferbot_to_top[1][1],NX*fdo2-1,MPI_FLOAT,INDEX[4],TAG4,MPI_COMM_WORLD,&req_send[3]);

	/* initialising of receive of buffer arrays. Same arrays for send and receive
	   are used. Thus, before starting receive, it is necessary to check if
	   Bsend has copied data into local buffers, i.e. has completed.
	*/
	MPI_Recv_init(&bufferlef_to_rig[1][1],NY*fdo2,MPI_FLOAT,INDEX[2],TAG1,MPI_COMM_WORLD,&req_rec[0]);
	MPI_Recv_init(&bufferrig_to_lef[1][1],NY*fdo2-1,MPI_FLOAT,INDEX[1],TAG2,MPI_COMM_WORLD,&req_rec[1]);
	MPI_Recv_init(&buffertop_to_bot[1][1],NX*fdo2,MPI_FLOAT,INDEX[4],TAG3,MPI_COMM_WORLD,&req_rec[2]);
	MPI_Recv_init(&bufferbot_to_top[1][1],NX*fdo2-1,MPI_FLOAT,INDEX[3],TAG4,MPI_COMM_WORLD,&req_rec[3]);

}
