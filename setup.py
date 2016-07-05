from setuptools import setup

setup(
	name = "motif_scraper",
	packages = ['motif_scraper'],
	version = "1.0.0",
	description = 'Tool for finding degenerate motifs in FASTA files',
	author = "Elisha Roberson",
	author_email = 'dr.eli.roberson@gmail.com',
	url = 'https://github.com/RobersonLab/motif_scraper',
	license = 'MIT',
	classifiers=[
	"Development Status :: 5 - Production/Stable",
	"Environment :: Console",
	"Intended Audience :: Science/Research",
	"License :: OSI Approved :: MIT License",
	"Topic :: Scientific/Engineering :: Bio-Informatics",
	"Programming Language :: Python :: 2.7",
	"Programming Language :: Python :: 3",
	"Programming Language :: Python :: 3.2",
	"Programming Language :: Python :: 3.3",
	"Programming Language :: Python :: 3.4",
	"Programming Language :: Python :: 3.5"
	],
	keywords='degenerate sequence motif site search',
	install_requires = ['pyfaidx', 'regex>=2016.01.10', 'six'],
	entry_points = {'console_scripts':["motif_scraper = motif_scraper.__main__:main"]},
	test_suite = 'nose.collector',
	tests_require = ['nose']
)
