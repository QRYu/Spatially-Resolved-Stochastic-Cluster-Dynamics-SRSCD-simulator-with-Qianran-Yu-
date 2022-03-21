// SRSCD -- Spatially Resolved Stochastic Cluster Dynamics
// Q Yu, January 2017.

#include<ctime>
#include<unistd.h>
#include<time.h>
#include"SCDWrapper.h"

int main() {
    //fork();
    SCDWrapper* srscd = new SCDWrapper(); /* establish spatially resolved scd */
    int64 theOtherKey = 0;
    Object* hostObject = nullptr;
    Reaction reaction = ERROR;
    int pointIndex = -1;
    int iStep = 0;
    double random;
    double advTime = 0.0;
    double dt = 0.0;
    double accTime = 0.0;
    long double bulkRate = 0.0; /* total rate of the whole bulk */
    int inputH = 0;
    double totalDPA = 0.2;
    double dpa = 0.0;
    fstream st;
    /* check whether to restart*/
    restart(iStep, advTime, srscd);
    srscd->getAndExamineRate();
    /* check ended */
    srand(time(0));
    /*display damage*/
    srscd->displayDamage();
    srscd->displayAllObject();
    srscd->drawSpeciesAndReactions(advTime);
    clock_t t0, t1;
    t0 = clock();
    //srscd->drawHD(advTime);
    //while (advTime < TOTAL_TIME) {
    while(dpa < totalDPA){
        hostObject = srscd->selectReaction(theOtherKey, reaction, pointIndex);/* choose an event */
        srscd->processEvent(reaction, hostObject, pointIndex, theOtherKey, advTime, accTime); /* process event */
        if(reaction == 6){
            accTime = 0.0;
        }
        
        if(iStep%PSTEPS == 0){
            t1 = clock()-t0;
            st.open("st.txt", ios::app);
            st << (float)t1/CLOCKS_PER_SEC << "  "<< dpa << endl;
            
            /*
            cout << "t = " << advTime << endl;
            cout <<"iStep = "<< iStep << endl;
            cout<<"dt= " << dt <<endl;
            cout<<"BulkRate = "<<bulkRate<<endl;
             */
            srscd->drawSpeciesAndReactions(advTime);
            srscd->drawDamage(advTime);
            srscd->writeFile(advTime, iStep);
            //srscd->writeVacancy();
            //srscd->drawHD(advTime);
            st.close();
        }
        bulkRate = srscd->getAndExamineRate(); /* calculate the bulk rate */
        ++iStep;
        do {
            random = (double)rand() / RAND_MAX;
        } while (random == 0);
        dt = (-1) / bulkRate*log(random);
        accTime +=dt;
        advTime += dt;
        dpa = srscd->getTotalDpa();
    }
    srscd->drawSpeciesAndReactions(advTime);
    srscd->drawDamage(advTime);
    srscd->writeVacancy();
    srscd->writeSinkFile(advTime, iStep);
    cout<<"dpa = "<<dpa<<endl;
    cout << "Finished, Bye" << endl;
    return 0;
}

