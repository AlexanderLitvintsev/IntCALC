object FormSelectDir: TFormSelectDir
  Left = 404
  Top = 152
  Width = 432
  Height = 186
  Caption = 'Выберите папку'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object DirectoryListBox1: TDirectoryListBox
    Left = 0
    Top = 0
    Width = 145
    Height = 129
    ItemHeight = 16
    TabOrder = 0
    OnChange = DirectoryListBox1Change
  end
  object DriveComboBox1: TDriveComboBox
    Left = 0
    Top = 136
    Width = 145
    Height = 19
    TabOrder = 1
    OnChange = DriveComboBox1Change
  end
  object Button1: TButton
    Left = 344
    Top = 8
    Width = 75
    Height = 25
    Caption = 'ОК'
    TabOrder = 2
    OnClick = Button1Click
  end
  object Button2: TButton
    Left = 344
    Top = 40
    Width = 75
    Height = 25
    Caption = 'Отмена'
    TabOrder = 3
    OnClick = Button2Click
  end
  object FileListBox1: TFileListBox
    Left = 152
    Top = 0
    Width = 185
    Height = 153
    ItemHeight = 13
    TabOrder = 4
  end
end
