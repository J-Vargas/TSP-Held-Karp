//============================================================================
// Name        : TSP_Held_Karp.cpp
// Author      : Yamile Vargas
// Description : Traveling Salesman Problem with Held and Karp Algorithm
//============================================================================

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stack>
#include <queue>
#include <stdlib.h>
using namespace std;

#define INFINIT 999999

FILE *fp;

typedef struct st_point {
	int  x;
	int  y;
} point;

typedef struct st_subset {
	int id;
	int num_bits;
} subset;

typedef struct st_cost {
	double cost_;
	int from_point;
	int from_subset;
} cost;

subset *subsetBinary;
int  *subsetIndexes;
point *points;
int count_points =0;
double total_length = 0;

/*WINDOW DISPLAY: starts*/
Display *display_ptr;
Screen *screen_ptr;
int screen_num;
char *display_name = NULL;
unsigned int display_width, display_height;
Window win;
int border_width, win_x, win_y;
unsigned int win_width, win_height;
XWMHints *wm_hints;
XClassHint *class_hints;
XSizeHints *size_hints;
XTextProperty win_name, icon_name;
char *win_name_string = "TSP with Held & Karp Algorithms - YPV";
char *icon_name_string = "Icon for Example Window";

XEvent report;

GC gc,  gc_red,gc_delete,  gc_gray, gc_tred, gc_tyellow, gc_tblue;
unsigned long valuemask = 0;
XGCValues gc_values,  gc_red_values, gc_blue_values, gc_gray_values,gc_white_values;
XGCValues gc_tyellow_values, gc_tred_values, gc_tblue_values;
Colormap color_map;
XColor tmp_color1, tmp_color2;

void setColors(Display* display_ptr, Window	win){
	gc = XCreateGC( display_ptr, win, valuemask, &gc_values);
	XSetForeground( display_ptr, gc, BlackPixel( display_ptr, screen_num ) );
	XSetLineAttributes( display_ptr, gc, 4, LineSolid, CapRound, JoinRound);

	/* and three other graphics contexts, to draw in yellow and red and grey white=WhitePixel(dis, screen);*/
	gc_delete = XCreateGC( display_ptr, win, valuemask, &gc_white_values);
	XSetLineAttributes(display_ptr, gc_delete, 2, LineSolid,CapRound, JoinRound);
	if( XAllocNamedColor( display_ptr, color_map, "white",
			&tmp_color1, &tmp_color2 ) == 0 )
	{printf("failed to get color white\n"); exit(-1);}
	else
		XSetForeground( display_ptr, gc_delete, WhitePixel( display_ptr, screen_num ) );

	gc_tyellow = XCreateGC( display_ptr, win, valuemask, &gc_tyellow_values);
	XSetLineAttributes(display_ptr, gc_tyellow, 2, LineSolid,CapRound, JoinRound);
	if( XAllocNamedColor( display_ptr, color_map, "yellow", &tmp_color1, &tmp_color2 ) == 0 )
	{printf("failed to get color yellow\n"); exit(-1);}
	else{
		XSetForeground( display_ptr, gc_tyellow, tmp_color1.pixel );
	}

	/* other graphics contexts red*/
	gc_red = XCreateGC( display_ptr, win, valuemask, &gc_red_values);
	XSetLineAttributes( display_ptr, gc_red, 1, LineSolid, CapRound, JoinRound);
	gc_tred = XCreateGC( display_ptr, win, valuemask, &gc_tred_values);
	XSetLineAttributes( display_ptr, gc_tred, 2, LineSolid, CapRound, JoinRound);
	if( XAllocNamedColor( display_ptr, color_map, "red", &tmp_color1, &tmp_color2 ) == 0 )
	{printf("failed to get color red\n"); exit(-1);}
	else
	{XSetForeground( display_ptr, gc_tred, tmp_color1.pixel );
	XSetForeground( display_ptr, gc_red, tmp_color1.pixel );
	}
	/* other graphics contexts red*/

	gc_gray = XCreateGC( display_ptr, win, valuemask, &gc_gray_values);
	XSetLineAttributes( display_ptr, gc_gray, 2, LineSolid, CapRound, JoinRound);
	if( XAllocNamedColor( display_ptr, color_map, "light grey", &tmp_color1, &tmp_color2 ) == 0 )
	{printf("failed to get color gray\n"); exit(-1);}
	else
		XSetForeground( display_ptr, gc_gray, tmp_color1.pixel );

	/*other graphics contexts grey*/
	gc_tblue = XCreateGC( display_ptr, win, valuemask, &gc_tblue_values);
	XSetLineAttributes( display_ptr, gc_tblue, 2, LineSolid, CapRound, JoinRound);
	if( XAllocNamedColor( display_ptr, color_map, "blue", &tmp_color1, &tmp_color2 ) == 0 )
	{printf("failed to get color blue\n"); exit(-1);}
	else{
		XSetForeground( display_ptr, gc_tblue, tmp_color1.pixel );
	}
}

void draw_path( stack <point> path, int col){
	point  temp; 		GC temp_col ;
	if(col==1)
		temp_col = gc_tblue;
	else if (col==0)
		temp_col = gc_tred;
	else if (col==2)
		temp_col = gc_tyellow;
	else if (col==3)
		temp_col = gc_gray;
	else
		temp_col = gc_tyellow;

	if(!path.empty()){

		temp = path.top();
		path.pop();
		while(!path.empty()){
			XDrawLine(display_ptr, win, temp_col, temp.x, temp.y, path.top().x,  path.top().y);
			temp = path.top();
			path.pop();
		}
	}
}
/*
 * The merge_sort algorithm sorts the array of subsets by the increasing number of ones (bites)
 * */
subset  * merge_sort(subset *node, int n){
	int i = 0;
	subset * left_list = NULL;
	subset * right_list = NULL;
	if(n == 1)
		return node;
	left_list = new subset[n/2];
	right_list = new subset[n-(n/2)];

	for (i = 0; i< n/2;i++){
		left_list[i]=node[i];
		right_list[i]=node[i+(n/2)];
	}
	if(n%2!=0)
		right_list[n/2]=node[n-1];

	/*Recursion*/
	right_list = merge_sort(right_list, (n - n/2));
	left_list = merge_sort(left_list, n/2);

	subset *final_list = new subset[n];

	i = 0; int j = 0; int fin = 0;
	while((i < n/2) and (j < (n - n/2))){
		if(left_list[i].num_bits  <= right_list[j].num_bits){
			final_list[fin] = left_list[i];
			fin++; i++;
		}else{
			final_list[fin] = right_list[j];
			fin++; j++;
		}
	}
	if(i < n/2){
		while(i < n/2){
			final_list[fin] = left_list[i];
			fin++; i++;
		}
	}else if(j < (n - n/2)){
		while(j < (n - n/2)){
			final_list[fin] = right_list[j];
			fin++; j++;
		}
	}
	left_list = 0;
	right_list = 0;
	return final_list;
}
point  * merge_sort_x(point *node, int n){
	int i = 0;
	point * left_list = NULL;
	point * right_list = NULL;
	if(n == 1)
		return node;
	left_list = new point[n/2];
	right_list = new point[n-(n/2)];

	for (i = 0; i< n/2;i++){
		left_list[i]=node[i];
		right_list[i]=node[i+(n/2)];
	}
	if(n%2!=0)
		right_list[n/2]=node[n-1];

	/*Recursion*/
	right_list = merge_sort_x(right_list, (n - n/2));
	left_list = merge_sort_x(left_list, n/2);

	point *final_list = new point[n];

	i = 0; int j = 0; int fin = 0;
	while((i < n/2) and (j < (n - n/2))){
		if(left_list[i].x  <= right_list[j].x){
			final_list[fin] = left_list[i];
			fin++; i++;
		}else{
			final_list[fin] = right_list[j];
			fin++; j++;
		}
	}
	if(i < n/2){
		while(i < n/2){
			final_list[fin] = left_list[i];
			fin++; i++;
		}
	}else if(j < (n - n/2)){
		while(j < (n - n/2)){
			final_list[fin] = right_list[j];
			fin++; j++;
		}
	}
	left_list = 0;
	right_list = 0;
	return final_list;
}
/*
 * Function count_num_bits: returns the number of ones (bites) that a integer contains
 */
unsigned int count_num_bits(unsigned int n)
{
	unsigned int count = 0;
	while(n)
	{
		count = (n & 1)==1? count + 1: count;
		n >>= 1;
	}
	return count;
}
/*	Function getSubsets: returns a sorted array with the subsets of
 *	number_points = 2^number_points*/
void getSubsets(int number_points){
	int i = 0;
	int number_points_=(int)pow(2,number_points-1);
	subsetBinary = new subset[number_points_];
	for (i=0; i < number_points_; i++){
		subsetBinary[i].id = i;
		subsetBinary[i].num_bits = count_num_bits(i);
	}
	//sort the numbers by the amount of ones they contain
	subsetBinary = merge_sort(subsetBinary, number_points_);
	subsetIndexes  = new int[number_points_];
	for(i=0;i < number_points_; i++){
		subsetIndexes[subsetBinary[i].id]=i;
	}
}

/* Function belongsTo: evaluates if a number a or b belongs to subset_i
 * */
bool belongsTo(int a,int b,subset subset_i){
	int representation_a = (int)pow(2,a-1);
	int representation_b = (int)pow(2,b-1);
	bool result = false;
	if(a==0)
		result = ((representation_b & subset_i.id) == representation_b);
	else
		result = ((representation_a & subset_i.id) == representation_a) or
		((representation_b & subset_i.id) == representation_b);
	return result;
}

/*Distance between two points*/
double distance_p(point p1, point p2){
	double l = 0;
	if(!(p1.x==p2.x and p1.y == p2.y))
		l = (sqrt( pow( p2.x - p1.x, 2 ) + pow(p2.y - p1.y, 2 )));
	return l;
}

/*  Function get_min_cost: returns the min cost of going
 *  from points[i] to points[j] through sub_sets[k]
 */
cost get_min_cost(int point_a, int point_b, subset subset_i, cost **c, point *pnt){
	cost min_distance, temp;
	min_distance.cost_ = temp.cost_ = INFINIT;
	int i=0;
	int b=0;//representation of the bit
	int temp_subset = 0;
	int bit_index = 0;

	/* if subset_i is the empty subset returns distance from point_a to point_b */
	if(subset_i.num_bits == 0){
		min_distance.cost_ = distance_p(pnt[point_a],pnt[point_b]);
		min_distance.from_point = point_a;
		min_distance.from_subset = subset_i.id;
	}
	/* returns distance from point_a to point_b going through to subset_i
	 * as long as point_a and point_b are not in subset_i.
	 * */
	else if (!belongsTo(point_a, point_b, subset_i) or (point_b == 0))
	{//evaluates if int point_a or int point_b are in subset_i
		if(subset_i.num_bits==1){
			b=log2(subset_i.id) + 1;
			min_distance.cost_ = distance_p(pnt[b], pnt[point_b])+ c[b-1][0].cost_ ;
			min_distance.from_point = b-1;
			min_distance.from_subset = 0;
		}else{
			bit_index = (int)pow(2,i);
			while( bit_index < subset_i.id){
				if((bit_index & subset_i.id ) == bit_index ){ // if the element bit_index belongs to the subset
					b=log2(bit_index) + 1;
					temp_subset = subsetIndexes[bit_index xor subset_i.id];
					temp.cost_ = distance_p(pnt[b], pnt[point_b])+ c[b-1][temp_subset].cost_ ;
					temp.from_point = b-1;
					temp.from_subset = temp_subset;
					if(temp.cost_ < min_distance.cost_)
						min_distance = temp;
				}
				i++;
				bit_index = (int)pow(2,i);
			}
		}
	}
	return min_distance;
}

stack <point>  update_path(stack <point> a_path,point edge_a, point edge_b){
	stack <point> updated_path, p;
	queue <point> tmp;
	p = a_path;
	point tmp1;
	tmp1 = p.top();
	p.pop();
	while(!p.empty()){
		if(((tmp1.x==edge_a.x and tmp1.y==edge_a.y) or (tmp1.x==edge_b.x and tmp1.y==edge_b.y))
				and ((p.top().x == edge_a.x and p.top().y == edge_a.y) or (p.top().x == edge_b.x and
						p.top().y == edge_b.y ))){
			while(!p.empty()){
				updated_path.push(p.top());
				p.pop();
			}
		}
		if(!p.empty()){
			tmp.push(tmp1); 	tmp1 = p.top(); p.pop();
		}
	}
	tmp.push(tmp1);
	updated_path.pop();
	while(!tmp.empty()){
		updated_path.push(tmp.front());
		tmp.pop();
	}
	return updated_path;
}
stack <point> find_draw_bridge(stack <point> path, stack <point> q , int temp_total_length){
	stack <point> updated_path, updated_temp_path;
	stack <point> temp_q;
	stack <point> p = path;

	point temp_p,temp_q_j,min_edge_p_to,min_edge_p_from,min_edge_q_to,min_edge_q_from;
	point edge_p_to,edge_p_from,edge_q_to,edge_q_from;
	double temp_min_distance = INFINIT; double min_distance = INFINIT;
	temp_p = p.top();
	p.pop();
	if(path.size() ==  1)
		p.push(temp_p);
	while(!p.empty()){
		temp_q = q;
		temp_q_j = temp_q.top();
		temp_q.pop();
		while(!temp_q.empty()){
			temp_min_distance = temp_total_length - distance_p(temp_p,p.top()) - distance_p(temp_q_j, temp_q.top());
			if((distance_p(temp_p,temp_q_j)+ distance_p(p.top(),temp_q.top())) <
					(distance_p(temp_p,temp_q.top()) + distance_p(p.top(),temp_q_j))){
				temp_min_distance = temp_min_distance + distance_p(temp_p,temp_q_j)+ distance_p(p.top(),temp_q.top());
				min_edge_p_from = temp_p; min_edge_p_to = temp_q_j;
				min_edge_q_from = p.top(); min_edge_q_to = temp_q.top();
			}else{
				temp_min_distance = temp_min_distance + distance_p(temp_p,temp_q.top()) + distance_p(p.top(),temp_q_j);
				min_edge_p_from = temp_p; min_edge_p_to = temp_q.top();
				min_edge_q_from = p.top(); min_edge_q_to = temp_q_j;
			}
			if(temp_min_distance < min_distance){
				min_distance = temp_min_distance; // the new total distance
				edge_p_from= min_edge_p_from; edge_p_to= min_edge_p_to;
				edge_q_from= min_edge_q_from; edge_q_to= min_edge_q_to;
			}
			temp_q_j = temp_q.top();
			temp_q.pop();
		}
		temp_p = p.top();
		p.pop();
	}
	/*draw new edges */
	XDrawLine(display_ptr, win, gc_red, edge_p_from.x, edge_p_from.y, edge_p_to.x,  edge_p_to.y);
	XDrawLine(display_ptr, win, gc_red, edge_q_from.x, edge_q_from.y, edge_q_to.x,  edge_q_to.y);
	/*Delete old edges*/
	if(path.size() == 1 or path.size() == 3)
		XDrawLine(display_ptr, win, gc_delete,edge_p_to.x,  edge_p_to.y, edge_q_to.x,  edge_q_to.y);
	else{
		XDrawLine(display_ptr, win, gc_delete,edge_p_from.x, edge_p_from.y, edge_q_from.x, edge_q_from.y);
		XDrawLine(display_ptr, win, gc_delete,edge_p_to.x,  edge_p_to.y, edge_q_to.x,  edge_q_to.y);
	}

	//Updates the path, by removing from it  the edge that was deleted in the bridge.
	if(path.size() == 16)
		updated_path = update_path(path,edge_p_from,edge_q_from);

	total_length = min_distance;

	return updated_path;
}

stack <point> tracePath(cost final_cost,cost **c, point *pt){
	stack <point> path;
	int point_i = 0;
	int subset_i = 0;
	path.push(pt[0]);
	path.push(pt[final_cost.from_point + 1]);
	point_i = final_cost.from_point;
	subset_i = final_cost.from_subset;
	int mytemp = 0;
	while(!(subset_i == 0)){
		mytemp = c[point_i][subset_i].from_point;
		path.push(pt[mytemp + 1]);
		subset_i = c[point_i][subset_i].from_subset;
		point_i = mytemp;
	}
	if(subset_i == 0)
		path.push(pt[0]);
	return path;
}

stack <point> held_karp(subset *sub_sets, point *point_, int numberPoints){
	int i, point_i,num_subsets, subset_i;
	stack <point> path;
	i = point_i = subset_i = 0;
	num_subsets =(int)pow(2,numberPoints-1);
	cost **costs;
	costs = new cost*[numberPoints-1];
	for(i=0; i<numberPoints-1 ;i++){
		costs[i]= new cost[num_subsets];
	}
	for(subset_i = 0; subset_i < num_subsets-1; subset_i++ ){
		for(point_i = 1; point_i < numberPoints; point_i++){
			costs[point_i-1][subset_i]= get_min_cost(0,point_i,sub_sets[subset_i],costs,point_);
		}
	}
	cost final_cost = get_min_cost(0,0,sub_sets[num_subsets-1],costs,point_);
	total_length = total_length + final_cost.cost_;

	return tracePath(final_cost,costs,point_);
}

void  compute_optimal_TSP(){
	int i = 0 ; int p = 0;
	stack <point> path, temp_path;
	if(count_points > 14){
		points = merge_sort_x(points,count_points);
		getSubsets(15);
		point *pt = new point[15];
		for(i=15; i <= (count_points - (count_points%15)); i = i+15){
			for(p=0;p<15;p++){
				pt[p]=points[i-p-1];
			}
			path = held_karp(subsetBinary,pt,15);
			if(!temp_path.empty())
				temp_path = find_draw_bridge(path, temp_path,total_length);
			else
				temp_path = path;
			draw_path(temp_path, (count_points-i)%4);
		}
		subsetBinary = 0;
		subsetIndexes = 0;
	}
	if(count_points%15 != 0 ){
		if(count_points%15 ==1 ){
			point point_a = points[count_points - 1];
			stack <point>  temp1;
			temp1.push(point_a);
			path = temp1;
		}else 	if(count_points%15 ==2 ){
			point point_a = points[count_points - 1];
			point point_b = points[count_points - 2];
			stack <point>  temp2;
			temp2.push(point_a);
			temp2.push(point_b);
			temp2.push(point_a);
			path = temp2;
			total_length = total_length + 2*distance_p(point_a, point_b);
		}else{
			getSubsets(count_points%15);
			point *pot = new point[count_points%15];
			int count = (count_points - (count_points%15));
			for(p=0;p<count_points%15;p++){
				pot[p]=points[count+p];
			}
			path = held_karp(subsetBinary,pot,count_points%15);
		}
		draw_path(path, 5);
		if(!temp_path.empty())
		{temp_path = find_draw_bridge(path, temp_path,total_length);}
	}
}


void read_file(FILE *t){
	int px, py;
	int min_x, min_y;
	min_x = 0;

	while(!feof(t)){
		fscanf(t, "%d %d \n", &px, &py);
		if(min_x == 0){
			win_width=win_x = px; win_height=win_y =py;
			min_x=1;
		}else{
			win_x = min(win_x, px);
			win_y = min(win_y, py);
			win_width = max((int)win_width, px);
			win_height = max((int)win_height ,py);
		}
		points[count_points].x = px;
		points[count_points].y = py;

		count_points++;
	}fclose(t);
	min_x = win_x - ((win_width - win_x )* 0.1);
	min_y = win_y - ((win_height - win_y)* 0.1);

	win_width = (win_width - win_x )* 1.1;
	win_height = (win_height - win_y)* 1.1;
	win_x = min_x;
	win_y = min_y;
}

int main(int argc, char *argv[])
{
	int right_click = 0; int j =0;
	points = new point[1000];

	if(argc > 1){
		fp = fopen(argv[1], "r");
		read_file(fp);
	}else
		printf("File was not given.\n");

	/* opening display: basic connection to X Server */
	if( (display_ptr = XOpenDisplay(display_name)) == NULL ){printf("Could not open display. \n"); exit(-1);}
	screen_num = DefaultScreen( display_ptr );
	screen_ptr = DefaultScreenOfDisplay( display_ptr );
	color_map  = XDefaultColormap( display_ptr, screen_num );
	display_width  = DisplayWidth( display_ptr, screen_num );
	display_height = DisplayHeight( display_ptr, screen_num );
	/* creating the window */
	border_width = 10;
	if(argc <2){  win_x = 100; win_y = 100; win_width = display_width/1.7; win_height = (int) (win_width /1.1);  }
	/*rectangular window*/

	win= XCreateSimpleWindow( display_ptr, RootWindow( display_ptr, screen_num),
			win_x, win_y, win_width, win_height, border_width,
			BlackPixel(display_ptr, screen_num),
			WhitePixel(display_ptr, screen_num) );

	/* now try to put it on screen, this needs cooperation of window manager */
	size_hints = XAllocSizeHints();
	wm_hints = XAllocWMHints();
	class_hints = XAllocClassHint();
	if( size_hints == NULL || wm_hints == NULL || class_hints == NULL )
	{ printf("Error allocating memory for hints. \n"); exit(-1);}

	size_hints -> flags = PPosition | PSize | PMinSize  ;
	size_hints -> min_width = 60;
	size_hints -> min_height = 60;

	XStringListToTextProperty( &win_name_string,1,&win_name);
	XStringListToTextProperty( &icon_name_string,1,&icon_name);

	wm_hints -> flags = StateHint | InputHint ;
	wm_hints -> initial_state = NormalState;
	wm_hints -> input = False;

	class_hints -> res_name = "x_use_example";
	class_hints -> res_class = "examples";

	XSetWMProperties( display_ptr, win, &win_name, &icon_name, argv, argc, size_hints, wm_hints, class_hints );

	/* what events do we want to receive */
	XSelectInput( display_ptr, win, ExposureMask | StructureNotifyMask | ButtonPressMask );

	/* finally: put window on screen */
	XMapWindow( display_ptr, win );
	XFlush(display_ptr);

	/* create graphics context, so that we may draw in this window */
	setColors(display_ptr,win);

	while(1)
	{
		XNextEvent( display_ptr, &report );
		switch( report.type )
		{
		case Expose:
			if(count_points > 0){
				for (j=0; j < count_points;  j++){
					XDrawRectangle(display_ptr, win, gc_tblue,points[j].x, points[j].y,2,2);
				}
			}
			break;
		case ConfigureNotify:
			/* This event happens when the user changes the size of the window*/
			win_width = report.xconfigure.width;
			win_height = report.xconfigure.height;
			break;
		case ButtonPress:
		{
			int x, y;
			x = report.xbutton.x;
			y = report.xbutton.y;
			if (report.xbutton.button == Button1){
				if(right_click ==0 ){
					points[count_points].x = x;
					points[count_points].y = y;
					count_points++;
					XDrawRectangle(display_ptr, win, gc_tblue,x, y,2,2);
				}
			}else{ //right click
				right_click++;
				if (right_click == 1){
					if( count_points > 0){
						cout << "Number of points =  "<< count_points << endl;
						compute_optimal_TSP();
						cout << "Total length  of  the tour = " << total_length <<endl;
					}else
						printf("There were zero points\n");
				}
				if(right_click > 1){
					XFreeGC(display_ptr, gc);
					XFreeGC(display_ptr, gc_gray);
					XFreeGC(display_ptr, gc_red);
					XDestroyWindow(display_ptr,win);
					XCloseDisplay(display_ptr);
				}
			}
		}
		break;
		default:
			/* this is a catch-all for other events; it does not do anything.
	              One could look at the report type to see what the event was */
			break;
		}
	}
	exit(0);
}
