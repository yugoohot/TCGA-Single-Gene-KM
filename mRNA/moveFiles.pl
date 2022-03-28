use strict;
use File::Copy;

my $newDir="files";
unless(-d $newDir)
{
	mkdir $newDir or die $!;
}

my @allFiles=glob("*");
foreach my $subDir(@allFiles)
{
	if((-d $subDir) && ($subDir ne $newDir))
	{
		opendir(SUB,"./$subDir") or die $!;
		while(my $file=readdir(SUB))
		{
			if($file=~/\.gz$/)
			{
				#`cp ./$subDir/$file ./$newDir`;
				copy("$subDir/$file","$newDir") or die "Copy failed: $!";
			}
		}
		close(SUB);
	}
}


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
