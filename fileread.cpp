#include "fileread.h"

void CDataCenter::process(const std::string filepath)
{
    std::ifstream istrm(filepath);
    if(!istrm.is_open())
    {
        return ;
    }

    std::ofstream ostrm(RESDATAPATH);

    if(!ostrm.is_open())
    {
        return;
    }

    char readbuffer[READBUFFSIZE];
    char *token = NULL;
    static int i =0;
    while(istrm.peek()!= EOF)
    {
        float64 ellipse[2][5];
        memset((float64 *)ellipse , 0, sizeof (float64 )*2*5);
#if(0)
        memset(readbuffer,0,READBUFFSIZE*sizeof (char));
        istrm.getline(readbuffer, READBUFFSIZE);
        char *data = STRTOK(readbuffer, " ",token);
         ellipse[0][0]= std::stod(data);//double a

        data = STRTOK(token, " ",token);
        ellipse[0][1]= std::stod(data);//double b

        data = STRTOK(token, " ",token);
        ellipse[0][2]= std::stod(data);//double theta

        data = STRTOK(token, " ",token);
        ellipse[0][3]= std::stod(data);//double centerx

        data = STRTOK(token, " ",token);
        ellipse[0][4]= std::stod(data);//double centery
        ellipse[0][1] = ellipse[0][0];

        memset(readbuffer,0,READBUFFSIZE*sizeof (char));
        istrm.getline(readbuffer, READBUFFSIZE);
         data = STRTOK(readbuffer, " ",token);
         ellipse[1][0]= std::stod(data);//double a

        data = STRTOK(token, " ",token);
        ellipse[1][1]= std::stod(data);//double b

        data = STRTOK(token, " ",token);
        ellipse[1][2]= std::stod(data);//double theta

        data = STRTOK(token, " ",token);
        ellipse[1][3]= std::stod(data);//double centerx

        data = STRTOK(token, " ",token);
        ellipse[1][4]= std::stod(data);//double centery

        ellipse[1][1] = ellipse[1][0];



        stEllipse ellipse1 = EllipseGenerating(  ellipse[0][3] ,ellipse[0][4] ,ellipse[0][0],ellipse[0][1],Deg2Rad(ellipse[0][2]));
        stEllipse ellipse2 = EllipseGenerating(  ellipse[1][3] ,ellipse[1][4] ,ellipse[1][0],ellipse[1][1],Deg2Rad(ellipse[1][2]));
        stIntersection solution =
                ellipse_intersection(ellipse[0][3] ,ellipse[0][4] ,ellipse[0][0],ellipse[0][1],Deg2Rad(ellipse[0][2]),
                                     ellipse[1][3] ,ellipse[1][4] ,ellipse[1][0],ellipse[1][1],Deg2Rad(ellipse[1][2]));

       ostrm<< ellipse1.homogeneous.A<<","
            << ellipse1.homogeneous.B<<","
            << ellipse1.homogeneous.C<<","
            << ellipse1.homogeneous.D<<","
            << ellipse1.homogeneous.E<<","
            << ellipse1.homogeneous.F<<","
            << ellipse2.homogeneous.A<<","
            << ellipse2.homogeneous.B<<","
            << ellipse2.homogeneous.C<<","
            << ellipse2.homogeneous.D<<","
            << ellipse2.homogeneous.E<<","
            << ellipse2.homogeneous.F<<","
            << solution.point[0].x<<","
            << solution.point[0].y<<","
            << solution.point[1].x<<","
            << solution.point[1].y<<","
            << solution.point[2].x<<","
            << solution.point[2].y<<","
            << solution.point[3].x<<","
            << solution.point[3].y<<","
            <<solution.realnum<<std::endl;
#endif
//        for(int k = 0 ;k< solution.realnum;k++)
//        {
//            if(!CheckIsPoint(&ellipse1, solution.point[k]) || !CheckIsPoint(&ellipse2, solution.point[k]))
//            {
//                printf("--> error index : %d",i);
//            }
//        }
//        if(solution.realnum == 0)
//        {
//           printf("-> %d \n ",i);
//        }

        i+=2;



    }
    istrm.close();
    ostrm.close();
}
