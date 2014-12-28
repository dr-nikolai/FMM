program TestFMM;

uses
  Forms,
  TestFrameWork,
  GUITestRunner,
  Test_uFMM in 'Test_uFMM.pas',
  uFMM in '..\uFMM.pas';

{$R *.RES}

begin
  Application.Initialize;
  GUITestRunner.RunRegisteredTests;
end.
