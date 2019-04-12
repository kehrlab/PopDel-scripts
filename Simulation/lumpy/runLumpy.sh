source activate lumpy
/usr/bin/time -f "1\t%e\t%U\t%S\t%M" -ao time/lumpy.time lumpy/0001.lumpy.sh
/usr/bin/time -f "10\t%e\t%U\t%S\t%M" -ao time/lumpy.time lumpy/0010.lumpy.sh
/usr/bin/time -f "20\t%e\t%U\t%S\t%M" -ao time/lumpy.time lumpy/0020.lumpy.sh
/usr/bin/time -f "30\t%e\t%U\t%S\t%M" -ao time/lumpy.time lumpy/0030.lumpy.sh
/usr/bin/time -f "40\t%e\t%U\t%S\t%M" -ao time/lumpy.time lumpy/0040.lumpy.sh
/usr/bin/time -f "50\t%e\t%U\t%S\t%M" -ao time/lumpy.time lumpy/0050.lumpy.sh
/usr/bin/time -f "60\t%e\t%U\t%S\t%M" -ao time/lumpy.time lumpy/0060.lumpy.sh
/usr/bin/time -f "70\t%e\t%U\t%S\t%M" -ao time/lumpy.time lumpy/0070.lumpy.sh
/usr/bin/time -f "80\t%e\t%U\t%S\t%M" -ao time/lumpy.time lumpy/0080.lumpy.sh
/usr/bin/time -f "90\t%e\t%U\t%S\t%M" -ao time/lumpy.time lumpy/0090.lumpy.sh
/usr/bin/time -f "100\t%e\t%U\t%S\t%M" -ao time/lumpy.time lumpy/0100.lumpy.sh
/usr/bin/time -f "200\t%e\t%U\t%S\t%M" -ao time/lumpy.time lumpy/0200.lumpy.sh
source deactivate
