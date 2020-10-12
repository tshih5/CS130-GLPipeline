#include "driver_state.h"
#include <cstring>
#include <climits>
#include <cfloat>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete[] image_color;
    delete[] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state &state, int width, int height)
{
    state.image_width = width;
    state.image_height = height;
    state.image_color = 0;
    state.image_color = new pixel[height * width];
    state.image_depth = new float[height * width];
    for (int i = 0; i < height * width; ++i)
    {
        state.image_color[i] = make_pixel(0, 0, 0);
        state.image_depth[i] = FLT_MAX;
    }
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state &state, render_type type)
{
    //std::cout<<"TODO: implement rendering."<<std::endl;
    switch (type)
    {
    case render_type::invalid:
    {
    }
    case render_type::indexed:
    {
        int vertex_ind = 0; //index in state.vertex_data that points to the first value in the vertex
        data_vertex vertex; //information about the vertex is here.
        const data_geometry *triangle[3];
        data_geometry temp_triangle[3];
        for (int i = 0; i < state.num_triangles; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                vertex.data = state.vertex_data + state.index_data[vertex_ind] * state.floats_per_vertex;
                temp_triangle[j].data = vertex.data;
                vertex_ind++;
                triangle[j] = &temp_triangle[j];
                state.vertex_shader(vertex, temp_triangle[j], state.uniform_data);
            }
            clip_triangle(state, triangle, 0);
        }
    }

    case render_type::triangle:
    {
        int vertex_ind = 0; //index in state.vertex_data that points to the first value in the vertex
        data_vertex vertex; //information about the vertex is here.
        const data_geometry *triangle[3];
        data_geometry temp_triangle[3];

        //X  Y  Z  R  G  B  X  Y  Z  R  G  B  X  Y  Z  R  G  B  X  Y  Z  R  G  B  X  Y  Z  R  G  B  X  Y  Z  R  G  B  X  Y  Z  R  G  B     vertex data
        //0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41    data index
        //1                 2                 3                 4                 5                 6                 7                    start of vertex
        //1                                                     2                                                     3					   start of triangle
        for (int i = 0; i <= state.num_vertices / 3; ++i)
        { //for each triangle
            for (int j = 0; j < 3; ++j)
            {                                                                      //for each vertex of the triangle,
                vertex.data = &state.vertex_data[vertex_ind];                      //point to the address of the start of the vertex
                temp_triangle[j].data = vertex.data;                               //add vertex to temp triangle vertex j
                triangle[j] = &temp_triangle[j];                                   //point triangle to address of temp triangle to get around const
                vertex_ind += state.floats_per_vertex;                             //increment vertex_ind by number of floats per vertex
                state.vertex_shader(vertex, temp_triangle[j], state.uniform_data); //call vertex shader
            }
            //rasterize_triangle(state, triangle);                                  //call rasterize triangle
            clip_triangle(state, triangle, 0); //clip triangle
        }
        break;
    }

    case render_type::fan:
    {
        data_vertex vertex;               //information about the vertex is here.
        const data_geometry *triangle[3]; //points to the current triangle that will be sent to clipping
        data_geometry temp_triangle[3];
        int vertex_ind = 0;
        int num_triangles = state.num_vertices - 2;
        for (int i = 0; i < num_triangles; ++i)
        { //for each triangle
            for (int j = 0; j < 3; ++j)
            { //for each vertex of the triangle,
                if (i != 0 && j == 0)
                { //after the first triangle, we do not recalculate the first vertex;
                    continue;
                }
                vertex.data = &state.vertex_data[vertex_ind];                      //point to the address of the start of the vertex
                temp_triangle[j].data = vertex.data;                               //add vertex to temp triangle vertex j
                triangle[j] = &temp_triangle[j];                                   //point triangle to address of temp triangle to get around const
                vertex_ind += state.floats_per_vertex;                             //increment vertex_ind by number of floats per vertex
                state.vertex_shader(vertex, temp_triangle[j], state.uniform_data); //call vertex shader
            }
            vertex_ind -= state.floats_per_vertex; //decrement vertex index to get new second vertex;
            clip_triangle(state, triangle, 0);     //clip triangle
        }
        break;
    }

    case render_type::strip:
    {
        data_vertex vertex;               //information about the vertex is here.
        const data_geometry *triangle[3]; //points to the current triangle that will be sent to clipping
        data_geometry temp_triangle[3];
        int vertex_ind = 0;
        int num_triangles = state.num_vertices - 2;
        for (int i = 0; i < num_triangles; ++i)
        { //for each triangle
            for (int j = 0; j < 3; ++j)
            {                                                                      //for each vertex of the triangle,
                vertex.data = &state.vertex_data[vertex_ind];                      //point to the address of the start of the vertex
                temp_triangle[j].data = vertex.data;                               //add vertex to temp triangle vertex j
                triangle[j] = &temp_triangle[j];                                   //point triangle to address of temp triangle to get around const
                vertex_ind += state.floats_per_vertex;                             //increment vertex_ind by number of floats per vertex
                state.vertex_shader(vertex, temp_triangle[j], state.uniform_data); //call vertex shader
            }
            vertex_ind -= 2 * state.floats_per_vertex; //decrement vertex index by 2 to get new first vertex
            clip_triangle(state, triangle, 0);         //clip triangle
        }
        break;
    }
    }
}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state &state, const data_geometry *in[3], int face)
{
    if (face == 6)
    {
        rasterize_triangle(state, in);
        return;
    }
    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    //clip against right, left, top, bottom, far, near
    //              0       1    2     3      4    5
    int sign = pow(-1, face);
    //axis we will be clipping against
    unsigned axis = face / 2;
    //stores whether each vertex is inside the clipping face;
    bool inside[3] = {false};
    //number of inside faces
    unsigned num_inside = 0;
    //calculate which vertices are inside as well as the number that are inside
    for (unsigned i = 0; i < 3; ++i)
    {
        inside[i] = sign * in[i]->gl_Position[axis] <= in[i]->gl_Position[3];
        if (inside[i])
        {
            num_inside++;
        }
    }
    //no vertices are in the clipping face
    if (num_inside == 0)
    {
        //if all vertices are outside draw nothing
        return;
    }
    else if (num_inside == 3)
    {
        //if all vertices are inside pass to next face
        clip_triangle(state, in, face + 1);
    }
    else if (num_inside == 1)
    {
        //if one vertex is inside, clip triangle to face, pass new vertices to next face
        //create a new triangle with two new vertices
        const data_geometry *triangle[3];
        data_geometry temp_triangle[3];
        data_vertex vertex;

        //get the vertex that is inside
        unsigned vertex_inside = 4;
        for (unsigned i = 0; i < 3; ++i)
        {
            if (inside[i])
            {
                vertex_inside = i;
                break;
            }
        }
        //calculate new triangle and clip it agains the next face
        for (unsigned i = 0; i < 3; ++i)
        {
            if (!inside[i])
            {
                //interpolate position of the vertex with the vertex that is inside
                //p1 = alpha * a + (1 - alpha) * b
                //	a					b
                //in[vertex_inside], in[i];
                float lambda = (sign * in[i]->gl_Position[3] - in[i]->gl_Position[axis]) / (in[vertex_inside]->gl_Position[axis] - (sign * in[vertex_inside]->gl_Position[3]) + (sign * in[i]->gl_Position[3]) - in[i]->gl_Position[axis]);
            }
            else
            {
                //for the inside vertex; keep it the same
                //temp_triangle[i] = in[i];
            }
        }
    }
    else
    {
        //two vertices are inside
        //clip_triangle(state,triangle1,face+1);
        //clip_triangle(state,triangle2,face+1);
    }

    clip_triangle(state, in, face + 1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state &state, const data_geometry *in[3])
{
    //std::cout<<"TODO: implement rasterization"<<std::endl;
    //          (width              *                       x                     ) / 2   + (width     -         1) / 2
    float a_i = (state.image_width * in[0]->gl_Position[0] / in[0]->gl_Position[3]) * 0.5 + (state.image_width - 1) * 0.5;
    float a_j = (state.image_height * in[0]->gl_Position[1] / in[0]->gl_Position[3]) * 0.5 + (state.image_height - 1) * 0.5;
    float b_i = (state.image_width * in[1]->gl_Position[0] / in[1]->gl_Position[3]) * 0.5 + (state.image_width - 1) * 0.5;
    float b_j = (state.image_height * in[1]->gl_Position[1] / in[1]->gl_Position[3]) * 0.5 + (state.image_height - 1) * 0.5;
    float c_i = (state.image_width * in[2]->gl_Position[0] / in[2]->gl_Position[3]) * 0.5 + (state.image_width - 1) * 0.5;
    float c_j = (state.image_height * in[2]->gl_Position[1] / in[2]->gl_Position[3]) * 0.5 + (state.image_height - 1) * 0.5;

    //calculate z coordinates
    float a_z = in[0]->gl_Position[2] / in[0]->gl_Position[3];
    float b_z = in[1]->gl_Position[2] / in[1]->gl_Position[3];
    float c_z = in[2]->gl_Position[2] / in[2]->gl_Position[3];

    //get the max and min height/width of the triangle to iterate through instead of the entire image

    int max_width = std::min(std::max(std::max(a_i, b_i), c_i), state.image_width - float(1.0)) + 1;
    int min_width = std::max(std::min(std::min(a_i, b_i), c_i), float(0.0)) + 1;

    int max_height = std::min(std::max(std::max(a_j, b_j), c_j), state.image_height - float(1.0)) + 1;
    int min_height = std::max(std::min(std::min(a_j, b_j), c_j), float(0.0)) + 1;

    //std::cout << max_width << " " << min_width << " " << max_height << " " << min_height << std::endl;
    //std::cout << state.image_width << " " << state.image_height << std::endl;

    float area = ((b_i * c_j - c_i * b_j) - (a_i * c_j - c_i * a_j) + (a_i * b_j - b_i * a_j));

    for (int i = min_width; i < max_width; ++i)
    {
        for (int j = min_height; j < max_height; ++j)
        {
            //initialize new data fragment and point it to array data
            float data[MAX_FLOATS_PER_VERTEX];
            data_fragment fragment{data};
            //for calculating depth for z buffering
            float depth;
            //calculate the barycentric coordinates of the points
            float alpha_p = ((b_i * c_j - c_i * b_j) - (i * c_j - c_i * j) + (i * b_j - b_i * j)) / area;
            float beta_p = ((i * c_j - c_i * j) - (a_i * c_j - c_i * a_j) + (a_i * j - i * a_j)) / area;
            float gamma_p = ((b_i * j - i * b_j) - (a_i * j - i * a_j) + (a_i * b_j - b_i * a_j)) / area;

            depth = alpha_p * a_z + beta_p * b_z + gamma_p * c_z;

            if (alpha_p >= 0.0 && beta_p >= 0.0 && gamma_p >= 0.0 && depth < state.image_depth[i + (j * state.image_width)])
            {
                state.image_depth[i + (j * state.image_width)] = depth;
                for (int k = 0; k < state.floats_per_vertex; ++k)
                {
                    switch (state.interp_rules[k])
                    {
                    case interp_type::flat:
                    {
                        //pixel recieves value stored at the first vertex of the triangle
                        //add fragment data
                        fragment.data[k] = in[0]->data[k];
                        break;
                    }
                    case interp_type::smooth:
                    {
                        //Vertex values are interpolated using perspective-correct interpolation.
                        //we have alpha' beta' and gamma', but now we need to calculate alpha, beta, gamma
                        //temporary barycentric coordinate
                        float alpha = alpha_p / in[0]->gl_Position[3];
                        float beta = beta_p / in[1]->gl_Position[3];
                        float gamma = gamma_p / in[2]->gl_Position[3];
                        //calculate k
                        float k_1 = alpha + beta + gamma;
                        //finish calculating coordinatees
                        alpha /= k_1;
                        beta /= k_1;
                        gamma /= k_1;
                        //add fragment data
                        fragment.data[k] = (alpha * in[0]->data[k]) + (beta * in[1]->data[k]) + (gamma * in[2]->data[k]);
                        break;
                    }
                    case interp_type::noperspective:
                    {
                        //Vertex values are interpolated using image-space barycentric coordinates.
                        //just use alpha_p, beta_p, and gamma_p and add fragment data;
                        fragment.data[k] = (alpha_p * in[0]->data[k]) + (beta_p * in[1]->data[k]) + (gamma_p * in[2]->data[k]);
                    }
                    case interp_type::invalid:
                    {
                        break;
                    }
                    }
                }
                //call fragment shader
                data_output output;
                state.fragment_shader(fragment, output, state.uniform_data);
                //set the pixel color to output of fragment shader
                state.image_color[i + (j * state.image_width)] = make_pixel(output.output_color[0] * 255, output.output_color[1] * 255, output.output_color[2] * 255);
            }
        }
    }
}
