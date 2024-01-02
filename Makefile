
build:
	CC=/usr/bin/musl-gcc  cargo build --release --target=x86_64-unknown-linux-musl
#cross build --target x86_64-unknown-linux-musl

push:
	rsync -avz target/x86_64-unknown-linux-musl/release/excord-lr   u249633@sug-login4.hgsc.bcm.edu:/stornext/snfs5/next-gen/scratch/zhengxc/workspace/stix/test_hg002/excord-lr

h:build
	target/x86_64-unknown-linux-musl/release/excord-lr --help

t:build 
	time target/x86_64-unknown-linux-musl/release/excord-lr -b test/mini_aln.sort.bam -o test/t 


tnomerge:build 
	time target/x86_64-unknown-linux-musl/release/excord-lr -n -b test/mini_aln.sort.bam -o test/t-nomerge 


te:
	time samtools view -b  test/mini_aln.sort.bam |  test/excord  --discordantdistance 500 /dev/stdin  >test/te

t2:build
	time target/x86_64-unknown-linux-musl/release/excord-lr -b test/full_aln.sort.bam -o test/tt 

tteclean:
	rm test/ttex

tte: 
	time samtools view -b  test/mini_aln.sort.bam |  \
	/mnt/e/projects/git_repo/excord/cmd/excord/excord \
	--discordantdistance 500 /dev/stdin  >test/ttex
ttes:
	cat test/ttex test/t | sort | grep -v 'align'  | grep -v 'align' | sort -k 1,1 -k 2,2 -k 5,5 -k 6,6 |tt |less -SN 


tte2: 
	time samtools view -b  test/mini_aln.sort.bam |  \
	test/excord \
	--discordantdistance 500 /dev/stdin  >test/ttex2
tte2s:
	cat test/ttex2 test/t | sort | grep -v 'align'  | grep -v 'align' | sort -k 1,1 -k 2,2 -k 5,5 -k 6,6 |tt |less -SN 

t2f:
	sudo /home/zhengxc/.cargo/bin/flamegraph  \
	-o my_flamegraph.svg  \
	-- \
	target/x86_64-unknown-linux-musl/release/excord-lr   \
	-b test/full_aln.sort.bam \
	-o test/tt 

t3:build
	time target/x86_64-unknown-linux-musl/release/excord-lr -d -Q 1 -b test/visibleData.bam -o test/t3 

t4:build
	time target/x86_64-unknown-linux-musl/release/excord-lr --verbose -Q 1 -b test/1.1865639.bam -o test/t4 


treadname:build
	time target/x86_64-unknown-linux-musl/release/excord-lr --verbose  -b test/c27f2be4-a026-47a6-af4c-46d2e91c1b1b.new.bam -o test/c27f2be4-a026-47a6-af4c-46d2e91c1b1b.new 


treadname2:build
	time target/x86_64-unknown-linux-musl/release/excord-lr --verbose -p 0 -b test/c27f2be4-a026-47a6-af4c-46d2e91c1b1b.new.bam -o test/c27f2be4-a026-47a6-af4c-46d2e91c1b1b.new 

cp:build
	scp -r -oHostKeyAlgorithms=+ssh-rsa target/x86_64-unknown-linux-musl/release/excord-lr    zhengxc@10.51.131.125:/workspace/zhengxc/project/stix/test/

tlargeInsCigar:build
	time target/x86_64-unknown-linux-musl/release/excord-lr --verbose  -b test/m54328_180922_235017_57738046_ccstest.bam  -o test/m54328_180922_235017_57738046_ccstest.bed


tmiss:build
	time target/x86_64-unknown-linux-musl/release/excord-lr --verbose  -b test/20-62184901-62185021.sam  -o test/20-62184901-62185021.bed

tmiss2:build
	time target/x86_64-unknown-linux-musl/release/excord-lr --verbose  -b test/20-62184901-62185021_large.bam  -o test/20-62184901-62185021_large.bed

	
thard:build
	time target/x86_64-unknown-linux-musl/release/excord-lr --verbose  -b test/hardclip.bam  -o test/hardclip.bed