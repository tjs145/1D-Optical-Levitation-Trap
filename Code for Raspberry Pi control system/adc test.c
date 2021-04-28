#include <stdio.h>
#include <bcm2835.h>
#include <sys/time.h>

#define SAMPLES 1000000

int get_raw(){

	uint8_t mosi[10] = { 0x01, 0x02, 0x00 };
	uint8_t miso[10] = { 0 };
	bcm2835_spi_transfernb(mosi, miso, 3);
	int raw = (miso[2] + (miso[1]%0x0f << 8));
	return raw;
}

int main(){
	if (!bcm2835_init())return 1;
	
	bcm2835_spi_begin();
	bcm2835_spi_setBitOrder(BCM2835_SPI_BIT_ORDER_MSBFIRST);      // default
	bcm2835_spi_setDataMode(BCM2835_SPI_MODE0);                   // default
	bcm2835_spi_setClockDivider(BCM2835_SPI_CLOCK_DIVIDER_64);    // ~ 4 MHz on Pi 2 (probably higher on Pi 4)
	bcm2835_spi_chipSelect(BCM2835_SPI_CS0);                      // default
	bcm2835_spi_setChipSelectPolarity(BCM2835_SPI_CS0, LOW);      // default
	
	struct timeval tv;
	struct timezone tz;
	
	gettimeofday(&tv, &tz);
	
	printf("Seconds: %d\n", tv.tv_sec);
	printf("Microseconds: %d\n", tv.tv_usec);
	double start = tv.tv_sec + 0.000001*tv.tv_usec;
	printf("start: %lf\n", start);
	
	FILE *data;
	FILE *times;
	data = fopen("adc_data.txt", "w+");
	times = fopen("adc_times.txt", "w+");
	
	
	for (int i=0; i<SAMPLES; i++){
		
	int raw = get_raw();
	char raw_str[10];
	sprintf(raw_str, "%d\n", raw);
	fputs(raw_str, data);
	
	gettimeofday(&tv, &tz);
	double t = tv.tv_sec + 0.000001*tv.tv_usec;
	char t_str[10];
	sprintf(t_str, "%d\n", t);
	fputs(t_str, times);
	
	
}

	fclose(data);
	fclose(times);
	
	gettimeofday(&tv, &tz);
	
	printf("Seconds: %d\n", tv.tv_sec);
	printf("Microseconds: %d\n", tv.tv_usec);
	double end = tv.tv_sec + 0.000001*tv.tv_usec;
	printf("End: %lf\n", end);
	
	double total = end - start;
	double rate = SAMPLES/total;
	printf("Sample rate: %lf\n", rate);
	
	
	bcm2835_spi_end();
	bcm2835_close();
	return 0;
}
