/* very basic timer */
#include <ctime>
#include <iostream>
#include <string>
class timer{
public:
  timer():diff_time_sec(0.),started(false),t1(std::time(0)){};
  void start(){t1 = std::time(0);};
  void stop(){
    if(started) {std::cout <<  "You forgot to start the timer dumbass, giving time after previous start or constructor"  << std::endl ;}
    t2 = std::time(0); 
    diff_time_sec = std::difftime(t2,t1);
    std::cout << " difference: "<< diff_time_sec << " seconds " << std::endl;
  };
  void stop(const std::string& message){
    std::cout<< message << " " << std::endl; 
    stop();
    started = false;
  }
private:
  time_t t1;
  time_t t2;  
  double diff_time_sec;
  bool started;
};
