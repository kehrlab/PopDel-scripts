source activate lumpy
/usr/bin/time -f "1\t%e\t%U\t%S\t%M" -ao /time/svtyper.time lumpy/0001.svtyper.sh
/usr/bin/time -f "10\t%e\t%U\t%S\t%M" -ao /time/svtyper.time lumpy/0010.svtyper.sh
/usr/bin/time -f "20\t%e\t%U\t%S\t%M" -ao /time/svtyper.time lumpy/0020.svtyper.sh
/usr/bin/time -f "30\t%e\t%U\t%S\t%M" -ao /time/svtyper.time lumpy/0030.svtyper.sh
/usr/bin/time -f "40\t%e\t%U\t%S\t%M" -ao /time/svtyper.time lumpy/0040.svtyper.sh
/usr/bin/time -f "50\t%e\t%U\t%S\t%M" -ao /time/svtyper.time lumpy/0050.svtyper.sh
/usr/bin/time -f "60\t%e\t%U\t%S\t%M" -ao /time/svtyper.time lumpy/0060.svtyper.sh
/usr/bin/time -f "70\t%e\t%U\t%S\t%M" -ao /time/svtyper.time lumpy/0070.svtyper.sh
/usr/bin/time -f "80\t%e\t%U\t%S\t%M" -ao /time/svtyper.time lumpy/0080.svtyper.sh
/usr/bin/time -f "90\t%e\t%U\t%S\t%M" -ao /time/svtyper.time lumpy/0090.svtyper.sh
/usr/bin/time -f "100\t%e\t%U\t%S\t%M" -ao /time/svtyper.time lumpy/0100.svtyper.sh
source deactivate
