#ifndef __PATCH_GRID_H_
#define __PATCH_GRID_H_

class PatchPerGrid {
	private:
		int ngrid; 
		uInt *grid_nsp;

	public:
		PatchPerGrid(int _ngrid){
			this->ngrid=_ngrid;
			this->grid_nsp = new uInt[_ngrid*_ngrid];
			for(int i = 0 ; i<_ngrid*_ngrid ; i++){
				this->grid_nsp[i]=0;
			}
		}

		~PatchPerGrid(){
			delete[] this->grid_nsp;
		}

		int n_cell(){
			return this->ngrid;
		}

		uInt get_grid(int i, int j){
			return this->grid_nsp[i + j*ngrid];
		}

		/*uInt get_grid_direct(int chosen){
			return this->grid_nsp[chosen];
		}*/

		void change_nsp(int i, int j, int nsp){
			this->grid_nsp[i+j*ngrid]= nsp;
		}

		void substract_one_nsp(int i, int j){
			this->grid_nsp[i+j*ngrid]-=1;
		}


		void add_one_sp( int i, int j){
			this->grid_nsp[i+j*ngrid]+=1;
		}

		vector<int> rand_choose_gridD(uInt totalpatch){
			double temptotal = totalpatch* get_random();
			double sum = 0.0;
			vector<int> result(2);
			for (int i = 0 ; i < ngrid ; i++){
				for (int j = 0 ; j < ngrid ; j++){
					assert(0<=this->grid_nsp[i+j*ngrid]);
					sum+= (double) this->grid_nsp[i+j*ngrid];
					if( temptotal <= sum ){
						result[0]= i;
						result[1]=j;
						goto end;
					}
				}
			}
			end:
			return result;
		}
};

PatchPerGrid *initializePatchGrid(Grid *g){
	int ncell= g->get_number_of_cells();
	uInt newsp =0;

	PatchPerGrid *ppgrid = new PatchPerGrid(ncell);
	
	for(int i= 0 ; i<ncell ; i++){
		for(int j=0 ; j<ncell ; j++){
			newsp = g->count_per_grid(i,j);
			ppgrid->change_nsp(i,j,newsp);
		}
	}
	return ppgrid;
}

void output_ppg(PatchPerGrid *ppg,Grid *g){
	for( int i=0; i<ppg->n_cell(); i++){
		for( int j=0; j<ppg->n_cell(); j++){
			cout<<" in grid " << i<< " and "<<j<<" there are "<< ppg->get_grid(i,j)<<" patches\n";
			cout << " if calculated directly from grid there are "<<g->count_per_grid(i,j)<<endl;
		}
	}
}

#endif // __PATCH_GRID_H_
